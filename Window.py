import tkinter


class Window:
    def __init__(self, main):
        self.main = main
        self.frame = tkinter.Frame(self.main, background="red")
        self.canvas = tkinter.Canvas(self.frame, width=600, height=600)
