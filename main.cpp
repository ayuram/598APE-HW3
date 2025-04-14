#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

#define G 6.6743
#define dt 0.001
#define Gdt (G * dt)

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

float tdiff(struct timeval *start, struct timeval *end) {
  return (end->tv_sec - start->tv_sec) + 1e-6 * (end->tv_usec - start->tv_usec);
}

unsigned long long seed = 100;

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

int nplanets;
int timesteps;
float THETA = 2.0;  // Opening angle for Barnes-Hut
#define EPS   1e-4  // Softening parameter to avoid divide-by-zero


typedef struct Planet {
   double mass;
   double x;
   double y;
   double vx;
   double vy;
} Planet;

// Each quadtree node covers a region [xMin,xMax] x [yMin,yMax].
typedef struct Node {
   // Bounding box of this node
   double xMin, xMax;
   double yMin, yMax;

   // Total mass and center-of-mass (for the entire region)
   double mass;
   double cmx, cmy; // center of mass x,y

   // If this node is a "leaf" containing exactly one planet, hasPlanet=1
   int hasPlanet;
   Planet p;

   // Children (quad split): NW, NE, SW, SE
   struct Node* nw;
   struct Node* ne;
   struct Node* sw;
   struct Node* se;
} Node;

/* Create an empty node representing the bounding box [xMin,xMax] x [yMin,yMax]. */
Node* createNode(double xMin, double xMax, double yMin, double yMax) {
    Node* node = (Node*)calloc(1, sizeof(Node));
    node->xMin = xMin;  node->xMax = xMax;
    node->yMin = yMin;  node->yMax = yMax;
    node->mass = 0.0;   // will accumulate as we insert
    node->cmx  = 0.0;
    node->cmy  = 0.0;
    // children are NULL initially
    node->nw = node->ne = node->sw = node->se = NULL;
    node->hasPlanet = 0;
    return node;
}

/* Subdivide a node into 4 child quadrants (nw, ne, sw, se). */
void subdivide(Node* node) {
    double midX = 0.5 * (node->xMin + node->xMax);
    double midY = 0.5 * (node->yMin + node->yMax);

    // Northwest
    node->nw = createNode(node->xMin, midX, midY, node->yMax);
    // Northeast
    node->ne = createNode(midX, node->xMax, midY, node->yMax);
    // Southwest
    node->sw = createNode(node->xMin, midX, node->yMin, midY);
    // Southeast
    node->se = createNode(midX, node->xMax, node->yMin, midY);
}

/* Insert a planet p into the quadtree rooted at 'node'. */
void insertPlanet(Node* node, Planet p) {
    // If this node doesn't yet contain a planet and has no children, store p here
    if (node->hasPlanet == 0 && node->nw == NULL) {
        node->hasPlanet = 1;
        node->p = p;
        node->mass = p.mass;
        node->cmx  = p.x;
        node->cmy  = p.y;
        return;
    }

    // Otherwise, if it already has a single planet but no children, we subdivide first
    if (node->nw == NULL && node->hasPlanet == 1) {
        subdivide(node);
        // Insert the planet that was here into the correct child
        Planet old = node->p;
        node->hasPlanet = 0;  // It's no longer a leaf
        node->p.mass = 0;     // Clear out

        insertPlanet(node, old);
    }

    // Now we know we have children. Insert 'p' into the correct quadrant
    double midX = 0.5 * (node->xMin + node->xMax);
    double midY = 0.5 * (node->yMin + node->yMax);

    if (p.x <= midX) {
        if (p.y <= midY) {
            // SW
            insertPlanet(node->sw, p);
        } else {
            // NW
            insertPlanet(node->nw, p);
        }
    } else {
        if (p.y <= midY) {
            // SE
            insertPlanet(node->se, p);
        } else {
            // NE
            insertPlanet(node->ne, p);
        }
    }

    // Update this node's mass & center-of-mass
    node->mass += p.mass;
    node->cmx  += p.mass * p.x;
    node->cmy  += p.mass * p.y;
}

/* A small BFS or DFS pass to finalize center-of-mass for internal nodes. */
void finalizeNode(Node* node) {
    if (!node) return;

    // If node->mass > 0, center of mass is (node->cmx / node->mass, node->cmy / node->mass).
    if (node->mass > 1e-12) {
        node->cmx /= node->mass;
        node->cmy /= node->mass;
    }

    // Recurse
    finalizeNode(node->nw);
    finalizeNode(node->ne);
    finalizeNode(node->sw);
    finalizeNode(node->se);
}

Node* buildQuadTree(Planet* planets, int n) {
    // 1) Find bounding box
    double xMin = planets[0].x, xMax = planets[0].x;
    double yMin = planets[0].y, yMax = planets[0].y;

    for (int i = 1; i < n; i++) {
        if (planets[i].x < xMin) xMin = planets[i].x;
        if (planets[i].x > xMax) xMax = planets[i].x;
        if (planets[i].y < yMin) yMin = planets[i].y;
        if (planets[i].y > yMax) yMax = planets[i].y;
    }

    // Expand bounding box a bit to avoid zero-size issues
    double dx = xMax - xMin;
    double dy = yMax - yMin;
    // If all planets are nearly the same point, expand artificially
    if (dx < 1e-10) dx = 1e-2;
    if (dy < 1e-10) dy = 1e-2;

    xMax = xMin + ((dx > dy) ? dx : dy); // keep it square
    yMax = yMin + ((dx > dy) ? dx : dy);

    // 2) Create root node
    Node* root = createNode(xMin, xMax, yMin, yMax);

    // 3) Insert each planet
    for (int i = 0; i < n; i++) {
        insertPlanet(root, planets[i]);
    }

    // 4) Finalize center of mass
    finalizeNode(root);

    return root;
}

/* Recursively compute the contribution of 'node' to the force on planet p. */
void addForce(Node* node, Planet* p, double* fx, double* fy) {
    if (!node || node->mass < 1e-12) return;

    double dx = node->cmx - p->x;
    double dy = node->cmy - p->y;
    double distSquared = (dx*dx + dy*dy + EPS);
    double invDist = fast_rsqrt(distSquared);

    // "Size" of the region (largest dimension)
    double regionSize = (node->xMax - node->xMin);
    // Barnes-Hut criterion
    if ((node->hasPlanet && node->mass > 0) || (regionSize * invDist < THETA)) {
        // Treat this entire node as a single body
        // Gravitational force: F = G * (m1*m2) / r^2
        // direction = (dx/dist, dy/dist)
        // => acceleration on p: a = F/m1 = G * m2 / r^2
        // => dvx = dt * a_x = dt * ( G * node->mass * dx / r^3 )
        double F = (G * p->mass * node->mass) * (invDist * invDist);
        // We only need the portion that goes into planet p's acceleration
        // so effectively: F / p->mass => G * node->mass / dist^2
        double invDist3 = invDist * invDist * invDist;
        double ax = G * node->mass * dx * invDist3;
        double ay = G * node->mass * dy * invDist3;

        // accumulate into fx, fy
        *fx += ax;
        *fy += ay;
    } else {
        // Recurse on children
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
    //  For each planet, compute net force from the tree
    //    Then update velocity & position
    #pragma omp parallel for
    for (int i = 0; i < nplanets; i++) {
        double fx = 0.0, fy = 0.0; // accumulators for acceleration
        addForce(root, &planets[i], &fx, &fy);

        // Update velocity
        planets[i].vx += dt * fx;
        planets[i].vy += dt * fy;

        // Update position
        planets[i].x  += dt * planets[i].vx;
        planets[i].y  += dt * planets[i].vy;
    }
}

int main(int argc, const char** argv){
   if (argc < 3) {
      printf("Usage: %s <nplanets> <timesteps>\n", argv[0]);
      return 1;
   }
   nplanets  = atoi(argv[1]);
   timesteps = atoi(argv[2]);
   // if a third argument is given, use it as THETA
    if (argc > 3) {
        THETA = atof(argv[3]);
    }

   // Allocate initial planet array
   Planet* planets = (Planet*)malloc(sizeof(Planet) * nplanets);

   // Random init
   for (int i=0; i<nplanets; i++) {
    planets[i].mass = randomDouble() * 10 + 0.2;
    planets[i].x    = ( randomDouble() - 0.5 ) * 100 * pow(1 + nplanets, 0.4);
    planets[i].y    = ( randomDouble() - 0.5 ) * 100 * pow(1 + nplanets, 0.4);
    planets[i].vx   = randomDouble() * 5 - 2.5;
    planets[i].vy   = randomDouble() * 5 - 2.5;
   }

   // Time it
   struct timeval start, end;
   gettimeofday(&start, NULL);
   Node* root = buildQuadTree(planets, nplanets);
   for (int i = 0; i < timesteps; i++) {
    nextBarnesHut(planets, root);
   }
   freeTree(root);

   gettimeofday(&end, NULL);
   printf("Total time %0.6f seconds, final location %f %f\n",
         tdiff(&start, &end), planets[nplanets-1].x, planets[nplanets-1].y);

   // Clean up
   free(planets);

   return 0;
}
