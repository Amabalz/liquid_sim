import pygame
import math
import random
# global variables
size = 64
dt = 1/120
iter = 4
scale = 10
maxsize = size**2 - 1

# from x, y to x in an array
def xyToX(x, y):
    if x + y*size > maxsize:
        return maxsize
    else:
        return x + y*size


# FLUID CLASS
class FluidMat():
    diff = 0
    visc = 0

    preDye = [0] * (size**2)
    dye = [0] * (size**2)

    velx = [0] * (size**2)
    vely = [0] * (size**2)

    preVelx = [0] * (size**2)
    preVely = [0] * (size**2)
    
    def __init__(self, diffusion, viscosity):
        self.diff = diffusion
        self.visc = viscosity

    def addDye(self, x, y, amount):
        self.dye[xyToX(x, y)] += amount

    def addVelocity(self, x, y, amountx, amounty):
        self.velx[xyToX(x, y)] += amountx
        self.vely[xyToX(x, y)] += amounty



# USED FUNCTIONS
# diffuse function
def diffuse (b, x, x0, diff):
    a = dt * diff * (size - 2) ** 2
    lin_solve(b, x, x0, a, 1 + 6 * a)

# set bounds function
def set_bnd(b, x):
    for i in range(1, size-1):
        if b == 2:
            x[xyToX(i, 0)] = -x[xyToX(i, 1)]
            x[xyToX(i, size-1)] = -x[xyToX(i, size-2)]
        else:
            x[xyToX(i, 0)] = x[xyToX(i, 1)]
            x[xyToX(i, size-1)] = x[xyToX(i, size-2)]

    for j in range(1, size-1):
        if b == 1:
            x[xyToX(0, j)] = -x[xyToX(1, j)]
            x[xyToX(size-1, j)] = -x[xyToX(size-2, j)]
        else:
            x[xyToX(0, j)] = x[xyToX(1  , j)]
            x[xyToX(size-1, j)] = x[xyToX(size-2, j)]

    
    x[xyToX(0, 0)] = 0.5 * (x[xyToX(1, 0)] + x[xyToX(0, 1)])
    x[xyToX(0, size-1)] = 0.5 * (x[xyToX(1, size-1)] + x[xyToX(0, size-2)])
    x[xyToX(size-1, 0)] = 0.5 * (x[xyToX(size-2, 0)] + x[xyToX(size-1, 1)])
    x[xyToX(size-1, size-1)] = 0.5 * (x[xyToX(size-2, size-1)] + x[xyToX(size-1, size-2)])

# linear solve function
def lin_solve(b, x, x0, a, c):
    cRecip = 1.0 / c;
    for j in range(1, size-1):
        for i in range(1, size-1):
            x[xyToX(i, j)] = (x0[xyToX(i, j)] + a*(x[xyToX(i+1, j)] + x[xyToX(i-1, j)] + x[xyToX(i, j+1)] + x[xyToX(i, j-1)])) * cRecip
    set_bnd(b, x)

# fix fluid flow function
def fixFlow(velx, vely, p, div):

    for j in range(1, size-1):
        for i in range(1, size-1):
            div[xyToX(i, j)] = -0.5 * (velx[xyToX(i+1, j)] - velx[xyToX(i-1, j)] + vely[xyToX(i, j+1)] - vely[xyToX(i, j-1)])/size;
            p[xyToX(i, j)] = 0;
    set_bnd(0, div); 
    set_bnd(0, p);
    lin_solve(0, p, div, 1, 6);
    for j in range(1, size-1):
        for i in range(1, size-1):
            velx[xyToX(i, j)] -= 0.5 * (p[xyToX(i+1, j)] - p[xyToX(i-1, j)]) * size
            vely[xyToX(i, j)] -= 0.5 * (p[xyToX(i, j+1)] - p[xyToX(i, j-1)]) * size
    set_bnd(1, velx)
    set_bnd(2, vely)

# advection function
def advect(b, d, d0, velx, vely):
    dtx = dt * (size - 2)
    dty = dt * (size - 2)
    Nfloat = size;
    
    for j in range(1, size-1):
        for i in range(1, size-1):
            jfloat = j
            ifloat = i
            tmp1 = dtx * velx[xyToX(i, j)]
            tmp2 = dty * vely[xyToX(i, j)]
            x    = ifloat - tmp1
            y    = jfloat - tmp2
            
            if x < 0.5:
                x = 0.5
            if x > Nfloat + 0.5:
                x = Nfloat + 0.5
            i0 = math.floor(x) 
            i1 = i0 + 1
            if y < 0.5:
                y = 0.5
            if y > Nfloat + 0.5:
                y = Nfloat + 0.5 
            j0 = math.floor(y)
            j1 = j0 + 1
        
            s1 = x - i0
            s0 = 1 - s1
            t1 = y - j0
            t0 = 1 - t1
            
            i0i = i0
            i1i = i1
            j0i = j0
            j1i = j1
            
            d[xyToX(i, j)] = s0 * (t0 * d0[xyToX(i0i, j0i)] + t1 * d0[xyToX(i0i, j1i)]) + s1 * (t0 * d0[xyToX(i1i, j0i)] + t1 * d0[xyToX(i1i, j1i)])
    set_bnd(b, d);

# time step function
def fluidStep(fluid):
    visc = fluid.visc
    diff = fluid.diff
    Vx = fluid.velx
    Vy = fluid.vely
    V0x = fluid.preVelx
    V0y = fluid.preVely
    preDye = fluid.preDye
    dye = fluid.dye
    
    diffuse(1, V0x, Vx, visc);
    diffuse(2, V0y, Vy, visc);
    
    fixFlow(V0x, V0y, Vx, Vy);
    
    advect(1, Vx, V0x, V0x, V0y);
    advect(2, Vy, V0y, V0x, V0y);
    
    fixFlow(Vx ,Vy, V0x, V0y);
    
    diffuse(0, preDye, dye, diff);
    advect(0, dye, preDye, Vx, Vy);

# render dye function
def renderDye(fluid):
    for i in range(0, size):
        for j in range(0, size):
            x = i * scale
            y = j * scale
            color = min(fluid.dye[xyToX(i,j)], 255)
            pygame.draw.rect(screen, (color, color, color), pygame.Rect(x, y, scale, scale))

# SIMULATION
# initialize a fluid
fluid = FluidMat(0, 0)

# pygame setup
pygame.init()
screen = pygame.display.set_mode((size*scale, size*scale))
clock = pygame.time.Clock()
running = True
dt = 1/120 #delta time 1/FPS
mousePrevPos = [0]*2

while running:
	# poll for events
    # pygame.QUIT event means the user clicked X to close your window
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False
        else:
            if event.type == pygame.MOUSEBUTTONDOWN or event.type == pygame.MOUSEMOTION :
                mouse_presses = pygame.mouse.get_pressed()
                if mouse_presses[0]:
                    pos = pygame.mouse.get_pos()
                    fluid.addDye(int(pos[0]/scale), int(pos[1]/scale), 100)
                    amountX = (pos[0] - mousePrevPos[0])*1.2
                    amountY = (pos[1] - mousePrevPos[1])*1.2
                    fluid.addVelocity(int(pos[0]/scale), int(pos[1]/scale), amountX, amountY)
                    mousePrevPos = pos

    # fill the screen with a color to wipe away anything from last frame
    screen.fill("black")

    # add a spring of fluid
    fluid.addDye(int(size/2), int(size/2), 125)
    fluid.addVelocity(int(size/2), int(size/2), random.randint(-50,50),random.randint(-50,50))

    # step through fluid function
    fluidStep(fluid)
    renderDye(fluid)

    # flip() the display to put your work on screen
    pygame.display.flip()

    # limits FPS to 120
    clock.tick(120)
pygame.quit()