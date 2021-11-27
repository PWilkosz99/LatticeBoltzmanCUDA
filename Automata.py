from Utilities import *


def generateMatrix(ImageMatrix, StateMatrix):
    for x in range(400):
        for y in range(600):  # fill with value only part of matrixes
            StateMatrix[x][y][0] = 1  # example
            if(randomByProbability(0.5)):  # 0-1
                ImageMatrix[y][x] = [255, 0, 0]
            else:  # TODO repair weird triangle inside random image
                ImageMatrix[x][y] = [0, 0, 0]
    return ImageMatrix, StateMatrix
