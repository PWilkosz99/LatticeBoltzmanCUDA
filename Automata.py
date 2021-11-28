from tkinter import Canvas
import numpy
from Utilities import *
from PIL import ImageTk, Image


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

def transitionRule(ImageMatrix, StateMatrix, canvas):
    newImageMatrix = numpy.zeros([600, 800, 3], dtype=numpy.uint8) 
    newStateMatrx = numpy.zeros([600, 800, 5], dtype=numpy.uint8) 
    for x in range(800): #TODO automatic range
        for y in range(600):
            newImageMatrix[y][x] = [0, 255, 0]
    image = ImageTk.PhotoImage(image=Image.fromarray(newImageMatrix))
    canvas.image=image
    canvas.create_image(0, 0, anchor="nw", image=image)
