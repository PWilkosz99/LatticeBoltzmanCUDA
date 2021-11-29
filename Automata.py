from tkinter import Canvas
import numpy
from Utilities import *
from PIL import ImageTk, Image


def generateMatrix(ImageMatrix, StateMatrix):

    # Generate starting particles
    for x in range(294):
        for y in range(600):  # fill with value only part of matrixes
            if(randomByProbability(0.05)):  # 0-1
                ImageMatrix[y][x] = [255, 0, 0]
                # on the start each generated patricle on each cell flow in only one dirtection
                StateMatrix[y][x][random.randint(0, 3)] = 1
            else:  # TODO repair weird triangle inside random image
                ImageMatrix[x][y] = [0, 0, 0]

    # Generate central wall
    for y in range(0, 250):
        for x in range(295, 305):
            ImageMatrix[y][x] = [255, 255, 255]
            StateMatrix[y][x][4] = 1
    for y in range(350, 600):
        for x in range(295, 305):
            ImageMatrix[y][x] = [255, 255, 255]
            StateMatrix[y][x][4] = 1
    return ImageMatrix, StateMatrix

# 0:N 1:S 2:W 3:E 4:SOLID


def transitionRule(ImageMatrix, StateMatrix, main, canvas):
    newImageMatrix = numpy.zeros(
        [ImageMatrix.shape[0], ImageMatrix.shape[1], 3], dtype=numpy.uint8)
    newStateMatrx = numpy.zeros(
        [ImageMatrix.shape[0], ImageMatrix.shape[1], 5], dtype=numpy.uint8)

    newImageMatrix, newStateMatrx = simulate(
        ImageMatrix, StateMatrix, newImageMatrix, newStateMatrx)

    image = ImageTk.PhotoImage(image=Image.fromarray(newImageMatrix))
    canvas.image = image
    canvas.create_image(0, 0, anchor="nw", image=image)
    main.after(1, lambda: transitionRule(newImageMatrix,
               newStateMatrx, main, canvas))  # recurensive


def simulate(ImageMatrix, StateMatrix, newImageMatrix, newStateMatrx):
    for x in range(newImageMatrix.shape[1]):  # X:1 Y:2
        for y in range(newImageMatrix.shape[0]):
            stateCouter = 0
            for s in range(5):
                if(StateMatrix[y][x][s] == 1):
                    stateCouter += 1
            if(stateCouter == 2 & StateMatrix[y][x][0] == 1 & StateMatrix[y][x][1] == 1 | stateCouter == 2 & StateMatrix[y][x][2] == 1 & StateMatrix[y][x][3] == 1):
                if(x <= 1 & y <= 1 & x >= newImageMatrix.shape[1]-1 & y >= newImageMatrix.shape[0]-1):
                    if(StateMatrix[y][x][0] == 1 & StateMatrix[y][x][1] == 1):
                        StateMatrix[y][x-1][2] = 1
                        StateMatrix[y][x+1][3] = 1
                        newImageMatrix[y][x-1] = [255, 0, 0]
                        newImageMatrix[y][x+2] = [255, 0, 0]
                    elif(StateMatrix[y][x][2] == 1 & StateMatrix[y][x][3] == 1):
                        StateMatrix[y-1][x][0] = 1
                        StateMatrix[y+1][x+1][1] = 1
                        newImageMatrix[y-1][x] = [255, 0, 0]
                        newImageMatrix[y+1][x] = [255, 0, 0]
                    #TODO : case when they bounce on 90'
            else:
                # wall stay on the same position and don't influence on other cells
                if(StateMatrix[y][x][4] == 1):
                    newStateMatrx[y][x][4] = 1
                    newImageMatrix[y][x] = [255, 255, 255]
                else:
                    if(StateMatrix[y][x][0] == 1):  # North
                        if(y <= 1 or StateMatrix[y-1][x][4] == 1):
                            newStateMatrx[y][x][1] = 1
                        else:
                            newStateMatrx[y-1][x][0] = 1
                            newImageMatrix[y-1][x] = [255, 0, 0]
                    if(StateMatrix[y][x][1] == 1):  # South
                        if(y >= newImageMatrix.shape[0]-1 or StateMatrix[y+1][x][4] == 1):
                            newStateMatrx[y][x][0] = 1
                        else:
                            newStateMatrx[y+1][x][1] = 1
                            newImageMatrix[y+1][x] = [255, 0, 0]
                    if(StateMatrix[y][x][2] == 1):  # West
                        if(x <= 1 or StateMatrix[y][x-1][2] == 1):
                            newStateMatrx[y][x][3] = 1
                        else:
                            newStateMatrx[y][x-1][2] = 1
                            newImageMatrix[y][x-1] = [255, 0, 0]
                    if(StateMatrix[y][x][3] == 1):  # East
                        if(x >= newImageMatrix.shape[1]-1 or StateMatrix[y][x+1][3] == 1):
                            newStateMatrx[y][x][2] = 1
                        else:
                            newStateMatrx[y][x+1][3] = 1
                            newImageMatrix[y][x+1] = [255, 0, 0]

    return newImageMatrix, newStateMatrx

# 0 1 2 3 *
# 1
# 2
# 3
# *
