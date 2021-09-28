import pygame as pg
from random import randint
from random import seed
import math

seed()

def hack_sig(x):
    if x >= 255:
        return 255
    if x <= 0:
        return 0
    return x
    # return 255 / (1 + math.exp(-x / 255))

k = 0.255
FPS = 200
BLOCK_SIZE = 1
B_WIDTH    = 512
B_HEIGHT   = 512

WIDTH=BLOCK_SIZE*B_WIDTH
HEIGHT=BLOCK_SIZE*B_WIDTH

class ColorDelta:
    def __init__(self, l):
        self._c = l.copy()

    def __iadd__(self, c):
        for i in range(3):
            self._c[i] += c._c[i]
        return self

    def __isub__(self, c):
        for i in range(3):
            self._c[i] -= c._c[i]
        return self

    def __imul__(self, a):
        for i in range(3):
            self._c[i] = a * self._c[i]
        return self

    def reset(self):
        for i in range(3):
            self._c[i] = 0

class Cell:
    def __init__(self, x, y):
        self._r = pg.Rect(x, y, BLOCK_SIZE, BLOCK_SIZE)
        r = randint(0,255)
        g = randint(0,255)
        b = randint(0,255)
        self._c = [r, g, b]
    
    def draw(self, sc):
        r = hack_sig(self._c[0])
        g = hack_sig(self._c[1])
        b = hack_sig(self._c[2])
        color = [r, g, b]
        
        pg.draw.rect(sc, color, self._r)

    def __sub__(self, c):
        return ColorDelta([self._c[i] - c._c[i] for i in range(3)])

    def __iadd__(self, c):
        for i in range(3):
            self._c[i] += c._c[i]
        return self

    def __isub__(self, c):
        for i in range(3):
            self._c[i] -= c._c[i]
        return self

class Grid:
    def __init__(self, alpha):
        self._g = [[Cell(j * BLOCK_SIZE, i * BLOCK_SIZE) for j in range(B_WIDTH)] for i in range(B_HEIGHT)]
        self._color_delta_matrix = [[ColorDelta([0.0,0.0,0.0]) for _ in range(B_WIDTH)] for _ in range(B_HEIGHT)]
        self._a = alpha

    def update(self):
        for i in range(B_HEIGHT):
            for j in range(B_WIDTH):
                d = self._g[i][(j+1) % B_WIDTH] - self._g[i][j]
                d *= self._a
                #
                self._color_delta_matrix[i][(j+1) % B_WIDTH] -= d
                self._color_delta_matrix[i][j]   += d
        #
        for j in range(B_WIDTH):
            for i in range(B_HEIGHT):
                d = self._g[(i+1) % B_HEIGHT][j] - self._g[i][j]
                d *= self._a
                #
                self._color_delta_matrix[(i+1) % B_HEIGHT][j] -= d
                self._color_delta_matrix[i][j]   += d
        #
        for i in range(B_HEIGHT):
            for j in range(B_WIDTH):
                self._g[i][j] += self._color_delta_matrix[i][j]
                self._color_delta_matrix[i][j].reset()

    def draw(self, sc):
        for i in range(B_HEIGHT):
            for j in range(B_WIDTH):
                self._g[i][j].draw(sc)

pg.init()
sc = pg.display.set_mode((WIDTH, HEIGHT))

running = True

clock = pg.time.Clock()

g = Grid(k)
g.draw(sc)
while running:
    #clock.tick(FPS)
    for event in pg.event.get():
        if event.type == pg.QUIT:
            running = False
    g.update()
    g.draw(sc)
    pg.display.flip()
