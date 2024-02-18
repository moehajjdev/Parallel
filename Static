#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <SDL2/SDL.h> 
#define MAX_ITER 1000
#define WIDTH 800
#define HEIGHT 800
typedef struct complex {
    double real;
    double imag;
} complex_t;
int mandelbrot(complex_t c) {
    complex_t z = {0, 0};
    int i;
    for (i = 0; i < MAX_ITER; i++) {
        if (sqrt(z.real * z.real + z.imag * z.imag) > 2) {
            return i;
        }
        complex_t tmp = {z.real * z.real - z.imag * z.imag + c.real, 2 * z.real * z.imag + c.imag};
        z = tmp;
    }
    return 0;
}
int main() {
    int i, j;
    static int data[WIDTH * HEIGHT];
    for (i = 0; i < WIDTH; i++) {
        for (j = 0; j < HEIGHT; j++) {
            complex_t c = {(double) i / WIDTH * 4 - 2, (double) j / HEIGHT * 4 - 2};
            int iter = mandelbrot(c);
            data[i * WIDTH + j] = iter;
        }
    }
    SDL_Window* window = SDL_CreateWindow("Mandelbrot Set", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, WIDTH, HEIGHT, 0);
    SDL_Surface* surface = SDL_GetWindowSurface(window);
    SDL_LockSurface(surface);
    for (i = 0; i < WIDTH; i++) {
        for (j = 0; j < HEIGHT; j++) {
            int index = i * WIDTH + j;
            Uint32 color = SDL_MapRGB(surface->format, 255 * (data[index] % 256) / 256, 255 * (data[index] % 256) / 256, 255 * (data[index] % 256) / 256);
            ((Uint32*) surface->pixels)[index] = color;
        }
    }
    SDL_UnlockSurface(surface);
    SDL_UpdateWindowSurface(window);
    SDL_Event event;
    while (SDL_WaitEvent(&event)) {
        if (event.type == SDL_QUIT) {
            break;
        }
    }
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
