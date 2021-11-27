import tkinter
import numpy
from Automata import *


class Window:
    def __init__(self, main):
        self.main = main
        self.ImageMatrix = numpy.zeros([600,800])
        self.OldMatrix = numpy.zeros([600,800])

        self.ImageMatrix, self.OldMatrix = generateMatrix(self.ImageMatrix, self.OldMatrix)
        self.frame = tkinter.Frame(self.main, background="red")
        self.canvas = tkinter.Canvas(self.frame, width=600, height=600)
