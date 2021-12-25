__constant__ double LATTICE_SPEED = 0.1;
__constant__ double TAU = 0.9;
__constant__ int DIRECTIONS = 9;
__constant__ int DIMENSIONS = 2;

#define KERNEL_HEADER(xvar, yvar, wvar, hvar) \
  const int x = blockIdx.x;\
  const int y = blockIdx.y;\
  const int width = gridDim.x;\
  const int height = gridDim.y
    

__global__ void streaming(double *out, double *in) {
  KERNEL_HEADER(x, y, width, height);

  int dir = 0;
  for(int i = -1; i <= 1; i++) {
    for(int j = -1; j <= 1; j++) {
        int target = y * DIRECTIONS + x * DIRECTIONS * height + dir;
        int xs = x + i;
        int ys = y + j;

        int source = ys * DIRECTIONS + xs * DIRECTIONS * height + dir;
        if(xs >= 0 && xs < width && ys >= 0 && ys < height) 
        {
            out[target] = in[source];
        } 
        else 
        {
            int bounce_dir = (1 - j) + (1 - i) * 3 ;
            int bounce_src = bounce_dir + y * DIRECTIONS + x * DIRECTIONS * height;
            out[target] = in[bounce_src];
        }
        dir++; 
    }
  }
}

__global__ void density(double *out, double* in) {
  KERNEL_HEADER(x, y, width, height);

  int target = y + x * height;
  out[target] = 0;

  int z = 0;
  for(int i = -1; i <= 1; i++) {
    for(int j = -1; j <= 1; j++) {
      int source = z + y * DIRECTIONS + x * DIRECTIONS * height;
      out[target] += in[source];
      z++;
    }
  }
}

__global__ void velocity(double *out, double* density, double* directions) {
  KERNEL_HEADER(x, y, width, height);

  // Compute target indices.
  const int target_x = 0 + y * DIMENSIONS + x * height * DIMENSIONS;
  const int target_y = 1 + y * DIMENSIONS + x * height * DIMENSIONS;
  out[target_x] = 0;
  out[target_y] = 0;

  int z = 0;
  for(int i = -1; i <= 1; i++) {
    for(int j = -1; j <= 1; j++) {
      const int source = z + y * DIRECTIONS + x * DIRECTIONS * height;
      out[target_x] += directions[source] * i;
      out[target_y] += directions[source] * j;
      z++;
    }
  }

  const int target = y + x * height;
  out[target_x] *= LATTICE_SPEED / density[target];
  out[target_y] *= LATTICE_SPEED / density[target];
}

__global__ void equilibrium(double *eq, double* density, double* velocity) {
  KERNEL_HEADER(x, y, width, height);

  const int density_src = y + x * height;
  const int xvel_src = 0 + y * DIMENSIONS + x * height * DIMENSIONS;
  const int yvel_src = 1 + y * DIMENSIONS + x * height * DIMENSIONS;

  // u . u
  double velmag = velocity[xvel_src] * velocity[xvel_src] + velocity[yvel_src] * velocity[yvel_src];

  int z = 0;
  for(int i = -1; i <= 1; i++) {
    for(int j = -1; j <= 1; j++) {
        // Compute the weight.
        double weight;
        if(i == 0 && j == 0) {
            weight = 4.0 / 9.0;
        } else if(i == 0 || j == 0) {
            weight = 1.0 / 9.0;
        } else {
            weight = 1.0 / 36.0;
        }

        // e_i . u
        double dotprod = i * velocity[xvel_src] + j * velocity[yvel_src];

        double sum = 1.0;
        sum += 3 / LATTICE_SPEED * dotprod;
        sum += 4.5 / (LATTICE_SPEED * LATTICE_SPEED) * dotprod * dotprod;
        sum -= 1.5 / (LATTICE_SPEED * LATTICE_SPEED) * velmag;

        const int target = z + y * DIRECTIONS + x * DIRECTIONS * height;
        eq[target] = weight * density[density_src] * sum;

        z++;
    }
  }
}

__global__ void update(double *out, double* equilibrium, double* directions) {
  KERNEL_HEADER(x, y, width, height);

  int z = 0;
  for(int i = -1; i <= 1; i++) {
    for(int j = -1; j <= 1; j++) {
        const int target = z + y * DIRECTIONS + x * DIRECTIONS * height;
        out[target] = directions[target] - (directions[target] - equilibrium[target]) / TAU;
        z++;
    }
  }
}
