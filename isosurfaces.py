#!/usr/bin/python

import vtk
import sys
import math
import numpy as np
from vtk.util import numpy_support as npvtk

#vtkarray = npvtk.numpy_to_vtk(numpy_array)
#numpy_array = npvtk.vtk_to_numpy(vtkarray)

# Create the RenderWindow, Renderer and both Actors
def Sphere():
    # create the quadric function definition
    quadric = vtk.vtkQuadric()
    quadric.SetCoefficients(1, 1, 1, 0, 0, 0, 1, 0, 0, 0)

    # F(x,y,z) = a0*x^2 + a1*y^2 + a2*z^2 + a3*x*y + a4*y*z + a5*x*z + a6*x + a7*y + a8*z + a9
    # F(x,y,z) = 1*x^2 + 1*y^2 + 1*z^2

    PlotFunction(quadric, 50.)


def PlotFunction(quadric, value):
    # sample the quadric function
    sample = vtk.vtkSampleFunction()
    sample.SetSampleDimensions(500, 500, 500)
    sample.SetImplicitFunction(quadric)
    # double xmin = 0, xmax=1, ymin=0, ymax=1, zmin=0, zmax=1
    xmin = -100
    xmax = 110
    ymin = -100
    ymax = 100
    zmin = -100
    zmax = 100
    sample.SetModelBounds(xmin, xmax, ymin, ymax, zmin, zmax)

    # create the 0 isosurface
    contours = vtk.vtkContourFilter()
    contours.SetInputConnection(sample.GetOutputPort())
    contours.GenerateValues(1, value, value)

    # map the contours to graphical primitives
    contourMapper = vtk.vtkPolyDataMapper()
    contourMapper.SetInputConnection(contours.GetOutputPort())
    contourMapper.SetScalarRange(0.0, 1.2)

    # create an actor for the contours
    contourActor = vtk.vtkActor()
    contourActor.SetMapper(contourMapper)

    # -- create a box around the function to indicate the sampling volume --

    # create outline
    outline = vtk.vtkOutlineFilter()
    outline.SetInputConnection(sample.GetOutputPort())

    # map it to graphics primitives
    outlineMapper = vtk.vtkPolyDataMapper()
    outlineMapper.SetInputConnection(outline.GetOutputPort())

    # create an actor for it
    outlineActor = vtk.vtkActor()
    outlineActor.SetMapper(outlineMapper)
    outlineActor.GetProperty().SetColor(0, 0, 0)

    # setup the window
    ren1 = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren1)
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    # add the actors to the scene
    ren1.AddActor(contourActor)
    ren1.AddActor(outlineActor)
    ren1.SetBackground(1, 1, 1)  # Background color white

    # render and interact
    renWin.Render()
    iren.Start()

nx = 201    # nx
ny = 251    # ny
nz = 16     # nz

qc = np.fromfile('data/qc0035', dtype=np.float32).reshape((nz, ny, nx))
qcvtk = npvtk.numpy_to_vtk(qc)
PlotFunction(qcvtk, 50.)

renderer = vtk.vtkRenderer()
renderWindow = vtk.vtkRenderWindow()
renderWindow.AddRenderer(renderer)

interactor = vtk.vtkRenderWindowInteractor()
interactor.SetRenderWindow(renderWindow)

renderer.SetBackground(0.1, 0.2, 0.4)
renderWindow.SetSize(600, 600)

contourRep = vtk.vtkOrientedGlyphContourRepresentation()
contourRep.GetLinesProperty().SetColor(1, 0, 0)  # set color to red

contourWidget = vtk.vtkContourWidget()
contourWidget.SetInteractor(interactor)
contourWidget.SetRepresentation(contourRep)
contourWidget.On()

for arg in sys.argv:
    if "-Shift" == arg:
        contourWidget.GetEventTranslator().RemoveTranslation(
            vtk.vtkCommand.LeftButtonPressEvent)
        contourWidget.GetEventTranslator().SetTranslation(
            vtk.vtkCommand.LeftButtonPressEvent,
            vtk.vtkWidgetEvent.Translate)
    elif "-Scale" == arg:
        contourWidget.GetEventTranslator().RemoveTranslation(
            vtk.vtkCommand.LeftButtonPressEvent)
        contourWidget.GetEventTranslator().SetTranslation(
            vtk.vtkCommand.LeftButtonPressEvent,
            vtk.vtkWidgetEvent.Scale)


pd = vtk.vtkPolyData()

points = vtk.vtkPoints()
lines = vtk.vtkCellArray()

for i in range(0, 21):
    angle = 2.0 * math.pi * i / 20.0
    points.InsertPoint(i, 0.1 * math.cos(angle),
                       0.1 * math.sin(angle), 0.0)
    lines.InsertNextCell(i)

pd.SetPoints(points)
pd.SetLines(lines)

contourWidget.Initialize(pd, 1)
contourWidget.Render()
renderer.ResetCamera()
renderWindow.Render()

interactor.Initialize()
interactor.Start()

contourWidget.Off()

