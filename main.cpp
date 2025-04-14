#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>
#include <time.h>

#define G 6.6743
#define dt 0.001
#define Gdt (G * dt)
#define EPS 1e-4  // Softening parameter to avoid divide-by-zero

int nplanets;
int timesteps;
float THETA = 2.0;  // Opening angle for Barnes-Hut

unsigned long long seed = 100;

/* Fast inverse square root (single Newton-Raphson iteration). */
double fast_rsqrt(double number) {
    long i;
    double x2, y;
    const double threehalfs = 1.5;

    x2 = number * 0.5;
    y = number;
    i = *(long*)&y;
    i = 0x5fe6eb50c7b537a9 - (i >> 1);
    y = *(double*)&i;
    y = y * (threehalfs - (x2 * y * y));
    return y;
}

/* Time difference in seconds. */
float tdiff(struct timeval *start, struct timeval *end) {
    return (end->tv_sec - start->tv_sec) + 1e-6 * (end->tv_usec - start->tv_usec);
}

/* Simple random number generators. */
unsigned long long randomU64() {
    seed ^= (seed << 21);
    seed ^= (seed >> 35);
    seed ^= (seed << 4);
    return seed;
}

double randomDouble() {
    unsigned long long next = randomU64();
    next >>= (64 - 26);
    unsigned long long next2 = randomU64();
    next2 >>= (64 - 26);
    return ((next << 27) + next2) / (double)(1LL << 53);
}

/* Planet and Quadtree node definitions. */
typedef struct Planet {
    double mass;
    double x, y;
    double vx, vy;
} Planet;

typedef struct Node {
    double xMin, xMax, yMin, yMax;
    double mass;
    double cmx, cmy;
    int hasPlanet;
    Planet p;
    struct Node *nw, *ne, *sw, *se;
} Node;

/* Quadtree helpers. */
Node* createNode(double xMin, double xMax, double yMin, double yMax) {
    Node* node = (Node*)calloc(1, sizeof(Node));
    node->xMin = xMin; node->xMax = xMax;
    node->yMin = yMin; node->yMax = yMax;
    return node;
}

void subdivide(Node* node) {
    double midX = 0.5 * (node->xMin + node->xMax);
    double midY = 0.5 * (node->yMin + node->yMax);
    node->nw = createNode(node->xMin, midX, midY, node->yMax);
    node->ne = createNode(midX, node->xMax, midY, node->yMax);
    node->sw = createNode(node->xMin, midX, node->yMin, midY);
    node->se = createNode(midX, node->xMax, node->yMin, midY);
}

void insertPlanet(Node* node, Planet p) {
    if (node->hasPlanet == 0 && node->nw == NULL) {
        node->hasPlanet = 1;
        node->p = p;
        node->mass = p.mass;
        node->cmx = p.x;
        node->cmy = p.y;
        return;
    }

    if (node->nw == NULL && node->hasPlanet == 1) {
        subdivide(node);
        Planet old = node->p;
        node->hasPlanet = 0;
        node->p.mass = 0;
        insertPlanet(node, old);
    }

    double midX = 0.5 * (node->xMin + node->xMax);
    double midY = 0.5 * (node->yMin + node->yMax);

    if (p.x <= midX) {
        if (p.y <= midY) insertPlanet(node->sw, p);
        else insertPlanet(node->nw, p);
    } else {
        if (p.y <= midY) insertPlanet(node->se, p);
        else insertPlanet(node->ne, p);
    }

    node->mass += p.mass;
    node->cmx += p.mass * p.x;
    node->cmy += p.mass * p.y;
}

void finalizeNode(Node* node) {
    if (!node) return;

    if (node->mass > 1e-12) {
        node->cmx /= node->mass;
        node->cmy /= node->mass;
    }

    finalizeNode(node->nw);
    finalizeNode(node->ne);
    finalizeNode(node->sw);
    finalizeNode(node->se);
}

Node* buildQuadTree(Planet* planets, int n) {
    double xMin = planets[0].x, xMax = planets[0].x;
    double yMin = planets[0].y, yMax = planets[0].y;

    for (int i = 1; i < n; i++) {
        if (planets[i].x < xMin) xMin = planets[i].x;
        if (planets[i].x > xMax) xMax = planets[i].x;
        if (planets[i].y < yMin) yMin = planets[i].y;
        if (planets[i].y > yMax) yMax = planets[i].y;
    }

    double dx = xMax - xMin;
    double dy = yMax - yMin;
    if (dx < 1e-10) dx = 1e-2;
    if (dy < 1e-10) dy = 1e-2;

    xMax = xMin + ((dx > dy) ? dx : dy);
    yMax = yMin + ((dx > dy) ? dx : dy);

    Node* root = createNode(xMin, xMax, yMin, yMax);

    for (int i = 0; i < n; i++) {
        insertPlanet(root, planets[i]);
    }

    finalizeNode(root);
    return root;
}

/* Force computation. */
void addForce(Node* node, Planet* p, double* fx, double* fy) {
    if (!node || node->mass < 1e-12) return;

    double dx = node->cmx - p->x;
    double dy = node->cmy - p->y;
    double distSquared = (dx * dx + dy * dy + EPS);
    double invDist = fast_rsqrt(distSquared);
    double regionSize = (node->xMax - node->xMin);

    if ((node->hasPlanet && node->mass > 0) || (regionSize * invDist < THETA)) {
        double invDist3 = invDist * invDist * invDist;
        double ax = G * node->mass * dx * invDist3;
        double ay = G * node->mass * dy * invDist3;
        *fx += ax;
        *fy += ay;
    } else {
        addForce(node->nw, p, fx, fy);
        addForce(node->ne, p, fx, fy);
        addForce(node->sw, p, fx, fy);
        addForce(node->se, p, fx, fy);
    }
}

void freeTree(Node* node) {
    if (!node) return;
    freeTree(node->nw);
    freeTree(node->ne);
    freeTree(node->sw);
    freeTree(node->se);
    free(node);
}

void nextBarnesHut(Planet* planets, Node* root) {
    #pragma omp parallel for
    for (int i = 0; i < nplanets; i++) {
        double fx = 0.0, fy = 0.0;
        addForce(root, &planets[i], &fx, &fy);

        planets[i].vx += dt * fx;
        planets[i].vy += dt * fy;
        planets[i].x += dt * planets[i].vx;
        planets[i].y += dt * planets[i].vy;
    }
}

/* Main simulation driver. */
int main(int argc, const char** argv) {
    if (argc < 3) {
        printf("Usage: %s <nplanets> <timesteps>\n", argv[0]);
        return 1;
    }

    nplanets = atoi(argv[1]);
    timesteps = atoi(argv[2]);

    if (argc > 3) {
        THETA = atof(argv[3]);
    }

    int max_threads = omp_get_max_threads();
    printf("Max threads: %d\n", max_threads);

    Planet* planets = (Planet*)malloc(sizeof(Planet) * nplanets);

    for (int i = 0; i < nplanets; i++) {
        planets[i].mass = randomDouble() * 10 + 0.2;
        planets[i].x = (randomDouble() - 0.5) * 100 * pow(1 + nplanets, 0.4);
        planets[i].y = (randomDouble() - 0.5) * 100 * pow(1 + nplanets, 0.4);
        planets[i].vx = randomDouble() * 5 - 2.5;
        planets[i].vy = randomDouble() * 5 - 2.5;
    }

    srand(time(NULL));

    double rebuildThreshold = 0.1;
    double init_xMin, init_xMax, init_yMin, init_yMax;
    init_xMin = init_xMax = planets[0].x;
    init_yMin = init_yMax = planets[0].y;

    for (int i = 1; i < nplanets; i++) {
        if (planets[i].x < init_xMin) init_xMin = planets[i].x;
        if (planets[i].x > init_xMax) init_xMax = planets[i].x;
        if (planets[i].y < init_yMin) init_yMin = planets[i].y;
        if (planets[i].y > init_yMax) init_yMax = planets[i].y;
    }

    double dx = init_xMax - init_xMin;
    double dy = init_yMax - init_yMin;
    if (dx < 1e-10) dx = 1e-2;
    if (dy < 1e-10) dy = 1e-2;
    double baseSize = (dx > dy) ? dx : dy;
    init_xMax = init_xMin + baseSize;
    init_yMax = init_yMin + baseSize;

    struct timeval start, end;
    gettimeofday(&start, NULL);

    Node* root = buildQuadTree(planets, nplanets);
    int rebuildCount = 0;

    for (int t = 0; t < timesteps; t++) {
        nextBarnesHut(planets, root);

        if (t % 10 == 0) {
            int needsRebuild = 0;
            int sampleSize = nplanets / 20;
            if (sampleSize < 5) sampleSize = 5;

            for (int s = 0; s < sampleSize; s++) {
                int j = rand() % nplanets;
                if (planets[j].x < init_xMin - rebuildThreshold * baseSize ||
                    planets[j].x > init_xMax + rebuildThreshold * baseSize ||
                    planets[j].y < init_yMin - rebuildThreshold * baseSize ||
                    planets[j].y > init_yMax + rebuildThreshold * baseSize) {
                    needsRebuild = 1;
                    break;
                }
            }

            if (needsRebuild) {
                rebuildCount++;
                freeTree(root);
                root = buildQuadTree(planets, nplanets);

                init_xMin = init_xMax = planets[0].x;
                init_yMin = init_yMax = planets[0].y;
                for (int j = 1; j < nplanets; j++) {
                    if (planets[j].x < init_xMin) init_xMin = planets[j].x;
                    if (planets[j].x > init_xMax) init_xMax = planets[j].x;
                    if (planets[j].y < init_yMin) init_yMin = planets[j].y;
                    if (planets[j].y > init_yMax) init_yMax = planets[j].y;
                }

                dx = init_xMax - init_xMin;
                dy = init_yMax - init_yMin;
                if (dx < 1e-10) dx = 1e-2;
                if (dy < 1e-10) dy = 1e-2;
                baseSize = (dx > dy) ? dx : dy;
                init_xMax = init_xMin + baseSize;
                init_yMax = init_yMin + baseSize;
            }
        }
    }

    freeTree(root);
    gettimeofday(&end, NULL);

    printf("Total tree rebuilds: %d\n", rebuildCount);
    printf("Total time %0.6f seconds, final location %f %f\n",
        tdiff(&start, &end), planets[nplanets-1].x, planets[nplanets-1].y);

    free(planets);
    return 0;
}
