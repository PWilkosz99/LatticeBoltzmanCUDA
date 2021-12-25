import numpy
from PIL import ImageTk, Image
import numba
from numba import prange
from Utilities import *


@numba.jit(nopython=True, parallel=True)
def generateMatrix(ImageMatrix, StateMatrix):

    # Generate starting particles
    for x in prange(294):
        for y in prange(600):  # fill with value only part of matrixes
            if(random.random() < 0.5):  # 0-1
                ImageMatrix[y][x] = [255, 0, 0]
                # on the start each generated patricle on each cell flow in only one dirtection
                StateMatrix[y][x][random.randint(0, 3)] = 1
            else:
                ImageMatrix[x][y] = [0, 0, 0]

    # Generate central wall
    for y in prange(0, 250):
        for x in prange(295, 305):
            ImageMatrix[y][x] = [255, 255, 255]
            StateMatrix[y][x][4] = 1
    for y in prange(350, 600):
        for x in prange(295, 305):
            ImageMatrix[y][x] = [255, 255, 255]
            StateMatrix[y][x][4] = 1

    for y in prange(0, 600):
        ImageMatrix[y][1] = [255, 255, 255]
        StateMatrix[y][1][4] = 1
        ImageMatrix[y][799] = [255, 255, 255]
        StateMatrix[y][798][4] = 1
    for x in prange(0, 800):
        ImageMatrix[1][x] = [255, 255, 255]
        StateMatrix[1][x][4] = 1
        ImageMatrix[598][x] = [255, 255, 255]
        StateMatrix[598][x][4] = 1

    return ImageMatrix, StateMatrix

# 0:N 1:S 2:W 3:E 4:SOLID


def transitionRule(ImageMatrix, StateMatrix, main, canvas):
    newImageMatrix = numpy.zeros(
        [ImageMatrix.shape[0], ImageMatrix.shape[1], 3], dtype=numpy.uint8)
    newStateMatrx = numpy.zeros(
        [ImageMatrix.shape[0], ImageMatrix.shape[1], 5], dtype=numpy.uint8)

    newImageMatrix, newStateMatrx = simulate(
        StateMatrix, newImageMatrix, newStateMatrx)
    #print(simulate.parallel_diagnostics(level=4))

    image = ImageTk.PhotoImage(image=Image.fromarray(newImageMatrix))
    canvas.image = image
    canvas.create_image(0, 0, anchor="nw", image=image)
    main.after(1, lambda: transitionRule(newImageMatrix,
                                         newStateMatrx, main, canvas))  # recursively


@numba.jit(nopython=True, parallel=True)
def simulate(StateMatrix, newImageMatrix, newStateMatrx):
    for x in prange(1, newImageMatrix.shape[1]-1):  # X:1 Y:0
        for y in prange(1, newImageMatrix.shape[0]-1):
            stateTmp = StateMatrix[y][x]
            if(stateTmp[4] != 1):
                stateN = StateMatrix[y-1][x][1]
                stateS = StateMatrix[y+1][x][0]
                stateW = StateMatrix[y][x-1][3]
                stateE = StateMatrix[y][x+1][2]
                stateSolidN = StateMatrix[y-1][x][4]
                stateSolidS = StateMatrix[y+1][x][4]
                stateSolidW = StateMatrix[y][x-1][4]
                stateSolidE = StateMatrix[y][x+1][4]

                if(stateN == 1 and stateS == 1 and stateW == 0 and stateE == 0):
                    newStateMatrx[y][x][2] = 1
                    newStateMatrx[y][x][3] = 1
                    newImageMatrix[y][x] = [255, 0, 0]
                elif(stateN == 0 and stateS == 0 and stateW == 1 and stateE == 1):
                    newStateMatrx[y][x][0] = 1
                    newStateMatrx[y][x][1] = 1
                    newImageMatrix[y][x] = [255, 0, 0]
                else:
                    if(stateN == 1):
                        if(stateSolidS == 1):
                            newStateMatrx[y][x][0] = 1
                        else:
                            newStateMatrx[y][x][1] = 1
                        newImageMatrix[y][x] = [255, 0, 0]
                    if(stateS == 1):
                        if(stateSolidN == 1):
                            newStateMatrx[y][x][1] = 1
                        else:
                            newStateMatrx[y][x][0] = 1
                        newImageMatrix[y][x] = [255, 0, 0]
                    if(stateW == 1):
                        if(stateSolidE == 1):
                            newStateMatrx[y][x][2] = 1
                        else:
                            newStateMatrx[y][x][3] = 1
                        newImageMatrix[y][x] = [255, 0, 0]
                    if(stateE == 1):
                        if(stateSolidW == 1):
                            newStateMatrx[y][x][3] = 1
                        else:
                            newStateMatrx[y][x][2] = 1
                        newImageMatrix[y][x] = [255, 0, 0]
            else:  # wall stay on the same position and don't influence on other cells
                newStateMatrx[y][x][4] = 1
                newImageMatrix[y][x] = [255, 255, 255]
    return newImageMatrix, newStateMatrx

# 0 1 2 3 *
# 1
# 2
# 3
# *
