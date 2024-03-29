import tkinter
from Window import *
from Automata import *


def main():
    root = tkinter.Tk()
    root.withdraw()
    top = tkinter.Toplevel(root)
    top.protocol("WM_DELETE_WINDOW", root.destroy)

    appWindow = Window(top)
    top.title("Lattice Boltzman Automata by Piotr Wilkosz")
    top.mainloop()


if __name__ == '__main__':
    main()
