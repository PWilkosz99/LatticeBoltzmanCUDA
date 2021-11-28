from Utilities import *


def generateMatrix(ImageMatrix, StateMatrix):

    # Generate starting particles
    for x in range(294):
        for y in range(600):  # fill with value only part of matrixes
            StateMatrix[x][y][0] = 1  # example
            if(randomByProbability(0.5)):  # 0-1
                ImageMatrix[y][x] = [255, 0, 0]
            else:  # TODO repair weird triangle inside random image
                ImageMatrix[x][y] = [0, 0, 0]

    # Generate central wall
    for y in range(0, 250):
        for x in range(295, 305):
            ImageMatrix[y][x] = [255, 255, 255]
            StateMatrix[y][x] = 5
    for y in range(350, 600):
        for x in range(295, 305):
            ImageMatrix[y][x] = [255, 255, 255]
            StateMatrix[y][x] = 5
    return ImageMatrix, StateMatrix
