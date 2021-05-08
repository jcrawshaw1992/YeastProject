#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

import os 
from os import path
import math

# 0.2, 0.4, 0.6, 0.8,
Scalling = [ '0.4', '0.6', '0.8','1','1.2', '1.4', '1.6', '1.8', '2', '2.2', '2.4', '2.6', '2.8', '3']   

for S in  Scalling:

    Directory = '/Volumes/Hardrive/Projects/FSI/VaryingEdgeLength/RawData/HorizontalLength_'+S+'/'
    OutputDirectory = '/Volumes/Hardrive/Projects/FSI/VaryingEdgeLength/Rawdata/HorizontalLength_'+S+'/'
    
    if path.isdir(OutputDirectory)==0:
        os.mkdir(OutputDirectory)

   # load state
    LoadState('/Volumes/Hardrive/Projects/FSI/VaryingEdgeLength/RawData/HorizontalLength_'+S+'/State.pvsm', DataDirectory='/Volumes/Hardrive/Projects/FSI/VaryingEdgeLength/RawData/HorizontalLength_'+S,
        wholegeometryvelocity_FileName=['/Volumes/Hardrive/Projects/FSI/VaryingEdgeLength/FlowFiles/HorizontalLength_'+S+'/wholegeometry-velocity_0.vtu', '/Volumes/Hardrive/Projects/FSI/VaryingEdgeLength/FlowFiles/HorizontalLength_'+S+'/wholegeometry-velocity_'+S+'.vtu', '/Volumes/Hardrive/Projects/FSI/VaryingEdgeLength/FlowFiles/HorizontalLength_'+S+'/wholegeometry-velocity_2.vtu', '/Volumes/Hardrive/Projects/FSI/VaryingEdgeLength/FlowFiles/HorizontalLength_'+S+'/wholegeometry-velocity_3.vtu', '/Volumes/Hardrive/Projects/FSI/VaryingEdgeLength/FlowFiles/HorizontalLength_'+S+'/wholegeometry-velocity_4.vtu', '/Volumes/Hardrive/Projects/FSI/VaryingEdgeLength/FlowFiles/HorizontalLength_'+S+'/wholegeometry-velocity_5.vtu', '/Volumes/Hardrive/Projects/FSI/VaryingEdgeLength/FlowFiles/HorizontalLength_'+S+'/wholegeometry-velocity_6.vtu', '/Volumes/Hardrive/Projects/FSI/VaryingEdgeLength/FlowFiles/HorizontalLength_'+S+'/wholegeometry-velocity_7.vtu', '/Volumes/Hardrive/Projects/FSI/VaryingEdgeLength/FlowFiles/HorizontalLength_'+S+'/wholegeometry-velocity_8.vtu', '/Volumes/Hardrive/Projects/FSI/VaryingEdgeLength/FlowFiles/HorizontalLength_'+S+'/wholegeometry-velocity_9.vtu', '/Volumes/Hardrive/Projects/FSI/VaryingEdgeLength/FlowFiles/HorizontalLength_'+S+'/wholegeometry-velocity_10.vtu'])

    # find source
    cellSize4 = FindSource('CellSize4')

    # set active source
    SetActiveSource(cellSize4)

    # find view
    renderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')
    # uncomment following to set a specific view size
    # renderView1.ViewSize = [1552, 381]

    # set active view
    SetActiveView(renderView1)

    # get display properties
    cellSize4Display = GetDisplayProperties(cellSize4, view=renderView1)

    # set active source
    SetActiveSource(cellSize4)

    # show data in view
    cellSize4Display = Show(cellSize4, renderView1)

    # find view
    spreadSheetView2 = FindViewOrCreate('SpreadSheetView2', viewtype='SpreadSheetView')
    # uncomment following to set a specific view size
    # spreadSheetView2.ViewSize = [400, 400]

    # set active view
    SetActiveView(spreadSheetView2)

    # find source
    cellSize1 = FindSource('CellSize1')

    # set active source
    SetActiveSource(cellSize1)

    # show data in view
    cellSize4Display_1 = Show(cellSize4, spreadSheetView2)

    # trace defaults for the display properties.
    cellSize4Display_1.FieldAssociation = 'Cell Data'

    # set active source
    SetActiveSource(cellSize4)

    # save data
    SaveData(OutputDirectory+'Inlet.txt', proxy=cellSize4, Precision=40,
        FieldDelimiter='\t',
        FieldAssociation='Cells')

    # set active view
    SetActiveView(renderView1)

    # find source
    clip2 = FindSource('Clip2')

    # set active source
    SetActiveSource(clip2)

    # toggle 3D widget visibility (only when running from the GUI)
    Show3DWidgets(proxy=clip2.ClipType)

    # get display properties
    clip2Display = GetDisplayProperties(clip2, view=renderView1)

    # hide data in view
    Hide(clip2, renderView1)

    # set active source
    SetActiveSource(clip2)

    # show data in view
    clip2Display = Show(clip2, renderView1)

    # hide data in view
    Hide(clip2, renderView1)

    # set active source
    SetActiveSource(cellSize1)

    # get display properties
    cellSize1Display = GetDisplayProperties(cellSize1, view=renderView1)

    # get color transfer function/color map for 'velocity'
    velocityLUT = GetColorTransferFunction('velocity')

    # set active view
    SetActiveView(spreadSheetView2)

    # show data in view
    cellSize1Display_1 = Show(cellSize1, spreadSheetView2)

    # trace defaults for the display properties.
    cellSize1Display_1.FieldAssociation = 'Cell Data'

    # set active source
    SetActiveSource(cellSize1)

    # hide data in view
    Hide(cellSize1, spreadSheetView2)

    # show data in view
    cellSize1Display_1 = Show(cellSize1, spreadSheetView2)

    # show data in view
    cellSize1Display_1 = Show(cellSize1, spreadSheetView2)

    # set active view
    SetActiveView(renderView1)

    # find source
    cellSize2 = FindSource('CellSize2')

    # hide data in view
    Hide(cellSize2, renderView1)

    # find source
    slice3 = FindSource('Slice3')

    # hide data in view
    Hide(slice3, renderView1)

    # find source
    cellSize3 = FindSource('CellSize3')

    # hide data in view
    Hide(cellSize3, renderView1)

    # set active view
    SetActiveView(spreadSheetView2)

    # save data
    SaveData(OutputDirectory+'Lower.txt', proxy=cellSize1, Precision=40,
        FieldDelimiter='\t',
        WriteTimeSteps=1,
        FieldAssociation='Cells')

    # set active view
    SetActiveView(renderView1)

    # set active source
    SetActiveSource(cellSize2)

    # show data in view
    cellSize2Display = Show(cellSize2, renderView1)

    # trace defaults for the display properties.
    cellSize2Display.Representation = 'Surface'
    cellSize2Display.ColorArrayName = ['CELLS', 'velocity']
    cellSize2Display.LookupTable = velocityLUT
    cellSize2Display.LineWidth = 1.5
    cellSize2Display.OSPRayScaleArray = 'Area'
    cellSize2Display.OSPRayScaleFunction = 'PiecewiseFunction'
    cellSize2Display.SelectOrientationVectors = 'None'
    cellSize2Display.ScaleFactor = 7.80487789597828e-06
    cellSize2Display.SelectScaleArray = 'None'
    cellSize2Display.GlyphType = 'Arrow'
    cellSize2Display.GlyphTableIndexArray = 'None'
    cellSize2Display.DataAxesGrid = 'GridAxesRepresentation'
    cellSize2Display.PolarAxes = 'PolarAxesRepresentation'
    cellSize2Display.GaussianRadius = 3.90243894798914e-06
    cellSize2Display.SetScaleArray = [None, '']
    cellSize2Display.ScaleTransferFunction = 'PiecewiseFunction'
    cellSize2Display.OpacityArray = [None, '']
    cellSize2Display.OpacityTransferFunction = 'PiecewiseFunction'

    # show color bar/color legend
    cellSize2Display.SetScalarBarVisibility(renderView1, True)

    # find source
    slice2 = FindSource('Slice2')

    # set active source
    SetActiveSource(slice2)

    # toggle 3D widget visibility (only when running from the GUI)
    Show3DWidgets(proxy=slice2.SliceType)

    # get display properties
    slice2Display = GetDisplayProperties(slice2, view=renderView1)

    # set active source
    SetActiveSource(slice2)

    # show data in view
    slice2Display = Show(slice2, renderView1)

    # show color bar/color legend
    slice2Display.SetScalarBarVisibility(renderView1, True)

    # hide data in view
    Hide(cellSize1, renderView1)

    # hide data in view
    Hide(cellSize2, renderView1)

    # set active source
    SetActiveSource(cellSize2)

    # show data in view
    cellSize2Display = Show(cellSize2, renderView1)

    # show color bar/color legend
    cellSize2Display.SetScalarBarVisibility(renderView1, True)

    # set active view
    SetActiveView(spreadSheetView2)

    # set active source
    SetActiveSource(cellSize1)

    # show data in view
    slice3Display = Show(slice3, spreadSheetView2)

    # trace defaults for the display properties.
    slice3Display.FieldAssociation = 'Cell Data'

    # show data in view
    cellSize2Display_1 = Show(cellSize2, spreadSheetView2)

    # trace defaults for the display properties.
    cellSize2Display_1.FieldAssociation = 'Cell Data'

    # set active source
    SetActiveSource(slice2)

    # set active source
    SetActiveSource(slice2)

    # show data in view
    slice2Display_1 = Show(slice2, spreadSheetView2)

    # trace defaults for the display properties.
    slice2Display_1.FieldAssociation = 'Cell Data'

    # show data in view
    slice2Display_1 = Show(slice2, spreadSheetView2)

    # set active view
    SetActiveView(renderView1)

    # set active view
    SetActiveView(spreadSheetView2)

    # show data in view
    cellSize2Display_1 = Show(cellSize2, spreadSheetView2)

    # set active source
    SetActiveSource(cellSize2)

    # save data
    SaveData(OutputDirectory+'Upper.txt', proxy=cellSize2, Precision=40,
        FieldDelimiter='\t',
        WriteTimeSteps=1,
        FieldAssociation='Cells')

    # set active view
    SetActiveView(renderView1)

    # set active source
    SetActiveSource(slice3)

    # toggle 3D widget visibility (only when running from the GUI)
    Show3DWidgets(proxy=slice3.SliceType)

    # get display properties
    slice3Display_1 = GetDisplayProperties(slice3, view=renderView1)

    # set active source
    SetActiveSource(slice3)

    # show data in view
    slice3Display_1 = Show(slice3, renderView1)

    # hide data in view
    Hide(cellSize2, renderView1)

    # hide data in view
    Hide(slice2, renderView1)

    # hide data in view
    Hide(cellSize4, renderView1)

    # set active source
    SetActiveSource(cellSize3)

    # get display properties
    cellSize3Display = GetDisplayProperties(cellSize3, view=renderView1)

    # set active view
    SetActiveView(spreadSheetView2)

    # set active source
    SetActiveSource(cellSize2)

    # show data in view
    cellSize3Display_1 = Show(cellSize3, spreadSheetView2)

    # trace defaults for the display properties.
    cellSize3Display_1.FieldAssociation = 'Cell Data'

    # set active source
    SetActiveSource(cellSize3)

    # save data
    SaveData(OutputDirectory+'Middel.txt', proxy=cellSize3, Precision=40,
        FieldDelimiter='\t',
        WriteTimeSteps=1,
        FieldAssociation='Cells')

    # set active source
    SetActiveSource(cellSize3)

    # set active source
    SetActiveSource(cellSize4)

    # destroy cellSize2
    Delete(cellSize2)
    del cellSize2

    # destroy slice2
    Delete(slice2)
    del slice2

    # hide data in view
    Hide(cellSize3, spreadSheetView2)

    # show data in view
    slice3Display = Show(slice3, spreadSheetView2)

    # destroy cellSize3
    Delete(cellSize3)
    del cellSize3

    # hide data in view
    Hide(slice3, renderView1)

    # show data in view
    clip2Display = Show(clip2, renderView1)

    # destroy slice3
    Delete(slice3)
    del slice3

    # find source
    slice4 = FindSource('Slice4')

    # set active source
    SetActiveSource(slice4)

    # toggle 3D widget visibility (only when running from the GUI)
    Show3DWidgets(proxy=slice4.SliceType)

    # destroy cellSize4
    Delete(cellSize4)
    del cellSize4

    # find source
    clip4 = FindSource('Clip4')

    # set active source
    SetActiveSource(clip4)

    # toggle 3D widget visibility (only when running from the GUI)
    Show3DWidgets(proxy=clip4.ClipType)

    # destroy slice4
    Delete(slice4)
    del slice4

    # find source
    wholegeometryvelocity_ = FindSource('wholegeometry-velocity_*')

    # set active source
    SetActiveSource(wholegeometryvelocity_)

    # destroy clip4
    Delete(clip4)
    del clip4

    # destroy cellSize1
    Delete(cellSize1)
    del cellSize1

    # find source
    slice1 = FindSource('Slice1')

    # destroy slice1
    Delete(slice1)
    del slice1

    # hide data in view
    Hide(clip2, renderView1)

    # show data in view
    wholegeometryvelocity_Display = Show(wholegeometryvelocity_, renderView1)

    # trace defaults for the display properties.
    wholegeometryvelocity_Display.Representation = 'Surface'
    wholegeometryvelocity_Display.ColorArrayName = [None, '']
    wholegeometryvelocity_Display.OSPRayScaleArray = 'velocity'
    wholegeometryvelocity_Display.OSPRayScaleFunction = 'PiecewiseFunction'
    wholegeometryvelocity_Display.SelectOrientationVectors = 'None'
    wholegeometryvelocity_Display.ScaleFactor = 0.000137170735979453
    wholegeometryvelocity_Display.SelectScaleArray = 'None'
    wholegeometryvelocity_Display.GlyphType = 'Arrow'
    wholegeometryvelocity_Display.GlyphTableIndexArray = 'None'
    wholegeometryvelocity_Display.DataAxesGrid = 'GridAxesRepresentation'
    wholegeometryvelocity_Display.PolarAxes = 'PolarAxesRepresentation'
    wholegeometryvelocity_Display.ScalarOpacityUnitDistance = 1.0282837327226e-05
    wholegeometryvelocity_Display.GaussianRadius = 6.85853679897264e-05
    wholegeometryvelocity_Display.SetScaleArray = [None, '']
    wholegeometryvelocity_Display.ScaleTransferFunction = 'PiecewiseFunction'
    wholegeometryvelocity_Display.OpacityArray = [None, '']
    wholegeometryvelocity_Display.OpacityTransferFunction = 'PiecewiseFunction'

    # destroy clip2
    Delete(clip2)
    del clip2

    # destroy wholegeometryvelocity_
    Delete(wholegeometryvelocity_)
    del wholegeometryvelocity_

    # get animation scene
    animationScene1 = GetAnimationScene()

    # update animation scene based on data timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()

    #### saving camera placements for all active views

    # current camera placement for renderView1
    renderView1.CameraPosition = [0.001104865032704903, 0.00029202449231280304, 0.0015275709209871214]
    renderView1.CameraFocalPoint = [0.000430941348895431, 1.94175008800814e-05, 2.8030626702684107e-09]
    renderView1.CameraViewUp = [0.09688141781351473, 0.9715478254938379, -0.21612221926641065]
    renderView1.CameraParallelScale = 0.000437852070564865

    #### uncomment the following to render all views
    # RenderAllViews()
    # alternatively, if you want to write images, you can use SaveScreenshot(...).


      # # find view
    # spreadSheetView1 = FindViewOrCreate('SpreadSheetView1', viewtype='SpreadSheetView')
    # # uncomment following to set a specific view size
    # # spreadSheetView1.ViewSize = [400, 400]

    # # destroy spreadSheetView1
    # Delete(spreadSheetView1)
    # del spreadSheetView1

    # # find view
    # renderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')
    # # uncomment following to set a specific view size
    # # renderView1.ViewSize = [1552, 381]

    # # destroy renderView1
    # Delete(renderView1)
    # del renderView1

    # # get active view
    # spreadSheetView2 = GetActiveViewOrCreate('SpreadSheetView')
    # # uncomment following to set a specific view size
    # # spreadSheetView2.ViewSize = [400, 400]

    # # destroy spreadSheetView2
    # Delete(spreadSheetView2)
    # del spreadSheetView2

   