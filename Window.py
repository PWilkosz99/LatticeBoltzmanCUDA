import tkinter
from tkinter.constants import S
import numpy
from PIL import ImageTk, Image
from Automata import *


class Window:
    def __init__(self, main):
        self.main = main
        self.ImageMatrix = numpy.zeros([600, 800, 3], dtype=numpy.uint8)  # RGB
        self.StateMatrix = numpy.zeros([600, 800, 5], dtype=numpy.uint8)  # NSWE or solid

        self.ImageMatrix, self.StateMatrix = generateMatrix(
            self.ImageMatrix, self.StateMatrix)
        self.img = ImageTk.PhotoImage(image=Image.fromarray(self.ImageMatrix))

        self.frame = tkinter.Frame(self.main, background="white")
        self.canvas = tkinter.Canvas(self.frame, width=800, height=600)
        self.canvas.pack(side="left", fill="y")
        self.canvas.create_image(0, 0, anchor="nw", image=self.img)
        self.frame.pack()   

        self.main.after(50, lambda: transitionRule(self.ImageMatrix, self.StateMatrix, self.main, self.canvas))
