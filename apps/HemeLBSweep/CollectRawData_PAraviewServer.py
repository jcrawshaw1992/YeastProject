#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

import os 
from os import path
import math


Directory = '/data/vascrem/testoutput/HemeLBSweep/Length_Variation/FlowFiles/'
# load state
LoadState('/data/vascrem/testoutput/HemeLBSweep/Length_Variation/FlowFiles/LowerState.pvsm')

# print "A"

# A = 0.9325189791142847
# B = 1.8109551126806118
# C = 2.6272663534520184
# D = 3.475951713466752
# E = 4.334920411905742

# Scalling = ['0.4', '0.6', '0.8', '1.2', '1.4', '1.6', '1.8', '2', '2.2', '2.4', '2.6', '2.8', '3']   

# for S in  Scalling:
#     C = 2.6272663534520184+ A*(float(S)-1)
#     B = 1.8109551126806118+ A*(float(S)-1)

#     D = 3.475951713466752 + A*(float(S)-1)
#     D = D + (C-B)*(float(S)-1)
#     X =  D - math.sqrt(3)/2
#     XPosition =str(X)


#     E = 4.334920411905742 + A*(float(S)-1)
#     E = E + (C-B)*(float(S)-1)
#     MidPoint = (E + D) /2

#     FrontPoint = A*0.8*float(S)

#     Directory = '/Volumes/Hardrive/Projects/FSI/VaryingEdgeLength/FlowFiles/HorizontalLength_'+S+'/'
#     OutputDirectory = '/Volumes/Hardrive/Projects/FSI/VaryingEdgeLength/rawdata/HorizontalLength_'+S+'/'

#     if path.isdir(OutputDirectory)==0:
#         os.mkdir(OutputDirectory)


#     wholegeometryvelocity_ = XMLUnstructuredGridReader(FileName=[Directory +'wholegeometry-velocity_0.vtu', Directory +'wholegeometry-velocity_1.vtu', Directory +'wholegeometry-velocity_2.vtu', Directory +'wholegeometry-velocity_3.vtu', Directory +'wholegeometry-velocity_4.vtu', Directory +'wholegeometry-velocity_5.vtu', Directory +'wholegeometry-velocity_6.vtu', Directory +'wholegeometry-velocity_7.vtu', Directory +'wholegeometry-velocity_8.vtu', Directory +'wholegeometry-velocity_9.vtu', Directory +'wholegeometry-velocity_10.vtu'])
#     wholegeometryvelocity_.CellArrayStatus = ['velocity']

#     print "B"
#     # set active source
#     SetActiveSource(wholegeometryvelocity_)
#     print "C"
#     # get active view
#     spreadSheetView1 = GetActiveViewOrCreate('SpreadSheetView')
#     # uncomment following to set a specific view size
#     # spreadSheetView1.ViewSize = [400, 400]
#     print "D"
#     # show data in view
#     wholegeometryvelocity_Display = Show(wholegeometryvelocity_, spreadSheetView1)

#     # trace defaults for the display properties.
#     wholegeometryvelocity_Display.FieldAssociation = 'Cell Data'



#     print "Got to here"

#     # find source
#     clip1 = FindSource('Clip1')
#     print "E"
#     # set active source
#     SetActiveSource(clip1)

#     # toggle 3D widget visibility (only when running from the GUI)
#     Show3DWidgets(proxy=clip1.ClipType)
#     print "F"
#     # Properties modified on clip1
#     clip1.Input = wholegeometryvelocity_

#     # Properties modified on clip1.ClipType
#     # ####################clip1.ClipType.Input = wholegeometryvelocity_
#     print "F"
#     # find source
#     wholegeometryvelocity__1 = FindSource('wholegeometry-velocity_*')
#     print "G"   
#     # set active source
#     SetActiveSource(wholegeometryvelocity__1)

#     # # destroy wholegeometryvelocity__1
#     # Delete(wholegeometryvelocity_0vtu)
#     # del wholegeometryvelocity_0vtu

#     # set active source
#     SetActiveSource(wholegeometryvelocity_)

#     # set active source
#     SetActiveSource(clip1)

#     # show data in view
#     clip1Display = Show(clip1, spreadSheetView1)

#     # trace defaults for the display properties.
#     clip1Display.FieldAssociation = 'Cell Data'

#     # show data in view
#     clip1Display = Show(clip1, spreadSheetView1)

#     # set active source
#     SetActiveSource(wholegeometryvelocity_)

#     # show data in view
#     wholegeometryvelocity_Display = Show(wholegeometryvelocity_, spreadSheetView1)

#     # show data in view
#     wholegeometryvelocity_Display = Show(wholegeometryvelocity_, spreadSheetView1)

#     # find view
#     renderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')
#     # uncomment following to set a specific view size
#     # renderView1.ViewSize = [1552, 472]

#     # set active view
#     SetActiveView(renderView1)

#     # set active source
#     SetActiveSource(clip1)

#     # show data in view
#     clip1Display_1 = Show(clip1, renderView1)

#     # get color transfer function/color map for 'velocity'
#     velocityLUT = GetColorTransferFunction('velocity')

#     # get opacity transfer function/opacity map for 'velocity'
#     velocityPWF = GetOpacityTransferFunction('velocity')

#     # trace defaults for the display properties.
#     clip1Display_1.Representation = 'Surface'
#     clip1Display_1.ColorArrayName = ['CELLS', 'velocity']
#     clip1Display_1.LookupTable = velocityLUT
#     clip1Display_1.Opacity = 0.32
#     clip1Display_1.OSPRayScaleArray = 'velocity'
#     clip1Display_1.OSPRayScaleFunction = 'PiecewiseFunction'
#     clip1Display_1.SelectOrientationVectors = 'None'
#     clip1Display_1.ScaleFactor = 2.92918281047605e-05
#     clip1Display_1.SelectScaleArray = 'None'
#     clip1Display_1.GlyphType = 'Arrow'
#     clip1Display_1.GlyphTableIndexArray = 'None'
#     clip1Display_1.DataAxesGrid = 'GridAxesRepresentation'
#     clip1Display_1.PolarAxes = 'PolarAxesRepresentation'
#     clip1Display_1.ScalarOpacityFunction = velocityPWF
#     clip1Display_1.ScalarOpacityUnitDistance = 1.60562133934638e-05
#     clip1Display_1.GaussianRadius = 1.46459140523802e-05
#     clip1Display_1.SetScaleArray = [None, '']
#     clip1Display_1.ScaleTransferFunction = 'PiecewiseFunction'
#     clip1Display_1.OpacityArray = [None, '']
#     clip1Display_1.OpacityTransferFunction = 'PiecewiseFunction'

#     # show color bar/color legend
#     clip1Display_1.SetScalarBarVisibility(renderView1, True)

#     # set active source
#     SetActiveSource(wholegeometryvelocity_)

#     # show data in view
#     wholegeometryvelocity_Display_1 = Show(wholegeometryvelocity_, renderView1)

#     # trace defaults for the display properties.
#     wholegeometryvelocity_Display_1.Representation = 'Surface'
#     wholegeometryvelocity_Display_1.ColorArrayName = [None, '']
#     wholegeometryvelocity_Display_1.OSPRayScaleArray = 'velocity'
#     wholegeometryvelocity_Display_1.OSPRayScaleFunction = 'PiecewiseFunction'
#     wholegeometryvelocity_Display_1.SelectOrientationVectors = 'None'
#     wholegeometryvelocity_Display_1.ScaleFactor = 8.117073448374868e-05
#     wholegeometryvelocity_Display_1.SelectScaleArray = 'None'
#     wholegeometryvelocity_Display_1.GlyphType = 'Arrow'
#     wholegeometryvelocity_Display_1.GlyphTableIndexArray = 'None'
#     wholegeometryvelocity_Display_1.DataAxesGrid = 'GridAxesRepresentation'
#     wholegeometryvelocity_Display_1.PolarAxes = 'PolarAxesRepresentation'
#     wholegeometryvelocity_Display_1.ScalarOpacityUnitDistance = 7.893904493889662e-06
#     wholegeometryvelocity_Display_1.GaussianRadius = 4.058536724187434e-05
#     wholegeometryvelocity_Display_1.SetScaleArray = [None, '']
#     wholegeometryvelocity_Display_1.ScaleTransferFunction = 'PiecewiseFunction'
#     wholegeometryvelocity_Display_1.OpacityArray = [None, '']
#     wholegeometryvelocity_Display_1.OpacityTransferFunction = 'PiecewiseFunction'

#     # set active source
#     SetActiveSource(clip1)

#     print "Position of lower vessel "
#     # Properties modified on clip1.ClipType
#     clip1.ClipType.Position = [X, 0,0]
#     clip1.ClipType.Rotation = [0.1329530866338892, 0.052629105706948384, 0.7825191447420148]
#     clip1.ClipType.Scale = [0.18285273598001567, 0.18210095188736175, 2.2565686673233536]

#     # show data in view
#     wholegeometryvelocity_Display_1 = Show(wholegeometryvelocity_, renderView1)

#     # update the view to ensure updated data information
#     spreadSheetView1.Update()

#     # update the view to ensure updated data information
#     renderView1.Update()

#     # find view
#     spreadSheetView2 = FindViewOrCreate('SpreadSheetView2', viewtype='SpreadSheetView')
#     # uncomment following to set a specific view size
#     # spreadSheetView2.ViewSize = [400, 400]

#     # update the view to ensure updated data information
#     spreadSheetView2.Update()

#     # # Rescale transfer function
#     # velocityLUT.RescaleTransferFunction(6.74577295068e-10, 2.35521580988)

#     # # Rescale transfer function
#     # velocityPWF.RescaleTransferFunction(6.74577295068e-10, 2.35521580988)

#     # find source
#     slice1 = FindSource('Slice1')

#     # set active source
#     SetActiveSource(slice1)

#     # toggle 3D widget visibility (only when running from the GUI)
#     Show3DWidgets(proxy=slice1.SliceType)

#     # get display properties
#     slice1Display = GetDisplayProperties(slice1, view=renderView1)

#     # get color transfer function/color map for 'velocity'
#     velocityLUT_1 = GetColorTransferFunction('velocity')

#     print "Position of lower vessel slice"
#     # Properties modified on slice1.SliceType
#     slice1.SliceType.Origin = [X, 0,0]
#     slice1.SliceType.Normal = [0.6878333863612577, 0.7250742488202442, 0.03394946692605872]

#     # update the view to ensure updated data information
#     spreadSheetView1.Update()

#     # update the view to ensure updated data information
#     renderView1.Update()

#     # update the view to ensure updated data information
#     spreadSheetView2.Update()

#     # get color transfer function/color map for 'Area'
#     areaLUT = GetColorTransferFunction('Area')

#     # # Rescale transfer function
#     # areaLUT.RescaleTransferFunction(1.11244066298e-19, 2.00226119202e-11)

#     # get opacity transfer function/opacity map for 'Area'
#     areaPWF = GetOpacityTransferFunction('Area')

#     # # Rescale transfer function
#     # areaPWF.RescaleTransferFunction(1.11244066298e-19, 2.00226119202e-11)

#     # # find source
#     calculator1 = FindSource('Calculator1')

#     # get display properties
#     calculator1Display = GetDisplayProperties(calculator1, view=renderView1)

#     # get separate color transfer function/color map for 'Vel'
#     separate_calculator1Display_VelLUT = GetColorTransferFunction('Vel', calculator1Display, separate=True)

#     # Rescale transfer function
#     separate_calculator1Display_VelLUT.RescaleTransferFunction(-0.300152474862, 0.273336049423)

#     # get separate opacity transfer function/opacity map for 'Vel'
#     separate_calculator1Display_VelPWF = GetOpacityTransferFunction('Vel', calculator1Display, separate=True)

#     # Rescale transfer function
#     separate_calculator1Display_VelPWF.RescaleTransferFunction(-0.300152474862, 0.273336049423)

#     # set active source
#     SetActiveSource(slice1)

#     # show data in view
#     slice1Display = Show(slice1, renderView1)

#     # show color bar/color legend
#     slice1Display.SetScalarBarVisibility(renderView1, True)

#     # hide data in view
#     Hide(clip1, renderView1)

#     # hide data in view
#     Hide(wholegeometryvelocity_, renderView1)

#     print "Fix this bit"
#     # Properties modified on slice1.SliceType
#     slice1.SliceType.Origin = [X, 0,0]
#     slice1.SliceType.Normal = [0.560605867878971, 0.8257625525771993, 0.06194568314946608]

#     # update the view to ensure updated data information
#     spreadSheetView1.Update()

#     # update the view to ensure updated data information
#     renderView1.Update()

#     # update the view to ensure updated data information
#     spreadSheetView2.Update()

#     # Rescale transfer function
#     velocityLUT_1.RescaleTransferFunction(6.74577295068e-10, 2.35521580988)

#     # get opacity transfer function/opacity map for 'velocity'
#     velocityPWF_1 = GetOpacityTransferFunction('velocity')

#     # Rescale transfer function
#     velocityPWF_1.RescaleTransferFunction(6.74577295068e-10, 2.35521580988)

#     print "Fix this "
#     # Properties modified on slice1.SliceType
#     slice1.SliceType.Normal = [0.7, 0.7, 0]

#     # update the view to ensure updated data information
#     spreadSheetView1.Update()

#     # update the view to ensure updated data information
#     renderView1.Update()

#     # update the view to ensure updated data information
#     spreadSheetView2.Update()

#     print "Fix this "
#     # Properties modified on slice1.SliceType
#     # slice1.SliceType.Normal = [0.86602540378, 0.86602540378, 0]

#     # update the view to ensure updated data information
#     # spreadSheetView1.Update()

#     # update the view to ensure updated data information
#     # renderView1.Update()

#     # update the view to ensure updated data information
#     # spreadSheetView2.Update()

#     # change representation type
#     slice1Display.SetRepresentationType('Surface')

#     # toggle 3D widget visibility (only when running from the GUI)
#     Hide3DWidgets(proxy=slice1.SliceType)

#     # set active source
#     SetActiveSource(slice1)

#     # set active source
#     SetActiveSource(calculator1)

#     # Rescale transfer function
#     calculator1Display.ScaleTransferFunction.RescaleTransferFunction(-0.299904835477, -0.00386421500791)

#     # Rescale transfer function
#     calculator1Display.OpacityTransferFunction.RescaleTransferFunction(-0.299904835477, -0.00386421500791)

#     # set active view
#     SetActiveView(spreadSheetView2)

#     # find source
#     cellSize1 = FindSource('CellSize1')

#     # # set active source
#     # SetActiveSource(cellSize1)

  

#     print "Save the lower size"
#     # save data
#     SaveData(OutputDirectory+'LowerSize.txt', proxy=cellSize1, Precision=10,
#         FieldDelimiter='\t',
#         WriteTimeSteps=1,
#         FieldAssociation='Cells')

#     print "Have saved the first piece of data --  the lower size "

#     # set active source
#     SetActiveSource(calculator1)

#     # set active source
#     SetActiveSource(cellSize1)

#     # show data in view
#     calculator1Display_1 = Show(calculator1, spreadSheetView2)

#     # trace defaults for the display properties.
#     calculator1Display_1.FieldAssociation = 'Cell Data'

#     # set active source
#     SetActiveSource(calculator1)

#     print "Save the lower velocity"
#     # save data
#     SaveData(OutputDirectory+'LowerSize.txt', proxy=cellSize1, Precision=10,
#         FieldDelimiter='\t',
#         WriteTimeSteps=1,
#         FieldAssociation='Cells')


#     # set active source
#     SetActiveSource(calculator1)

#     # set active source
#     SetActiveSource(wholegeometryvelocity_)

#     # show data in view
#     wholegeometryvelocity_Display_2 = Show(wholegeometryvelocity_, spreadSheetView2)

#     # trace defaults for the display properties.
#     wholegeometryvelocity_Display_2.FieldAssociation = 'Cell Data'

#     # set active view
#     SetActiveView(renderView1)

#     # show data in view
#     wholegeometryvelocity_Display_1 = Show(wholegeometryvelocity_, renderView1)

#     # set active source
#     SetActiveSource(clip1)

#     # set active source
#     SetActiveSource(clip1)

#     # show data in view
#     clip1Display_1 = Show(clip1, renderView1)

#     # show color bar/color legend
#     clip1Display_1.SetScalarBarVisibility(renderView1, True)

#     # hide data in view
#     Hide(wholegeometryvelocity_, renderView1)

#     # set active source
#     SetActiveSource(slice1)

#     # set active source
#     SetActiveSource(clip1)

#     # Properties modified on clip1.ClipType
#     clip1.ClipType.Position = [X, 0,0]
#     clip1.ClipType.Rotation = [-1.0391275242407447, -0.20971833525893527, 1.5066717499058213]
#     clip1.ClipType.Scale = [0.1828527359800138, 0.18210095188736095, 2.256568667323315]

#     # update the view to ensure updated data information
#     spreadSheetView1.Update()
#     renderView1.Update()
#     spreadSheetView2.Update()

#     # set active source
#     SetActiveSource(slice1)

#     # toggle 3D widget visibility (only when running from the GUI)
#     Show3DWidgets(proxy=slice1.SliceType)


#     # Properties modified on slice1.SliceType
#     slice1.SliceType.Origin = [X, 0,0]
#     slice1.SliceType.Normal = [0.7062040601475004, -0.7062040601475004, 0.05051386802029261]

#     # update the view to ensure updated data information
#     spreadSheetView1.Update()
#     renderView1.Update()
#     spreadSheetView2.Update()

#     # Rescale transfer function
#     separate_calculator1Display_VelLUT.RescaleTransferFunction(-0.300152474862, 0.308601160326)

#     # Rescale transfer function
#     separate_calculator1Display_VelPWF.RescaleTransferFunction(-0.300152474862, 0.308601160326)

#     # set active source
#     SetActiveSource(wholegeometryvelocity_)

#     # show data in view
#     wholegeometryvelocity_Display_1 = Show(wholegeometryvelocity_, renderView1)

#     # set active source
#     SetActiveSource(slice1)

#     # hide data in view
#     Hide(clip1, renderView1)

#     # hide data in view
#     Hide(wholegeometryvelocity_, renderView1)

#     # set active source
#     SetActiveSource(wholegeometryvelocity_)

#     # show data in view
#     wholegeometryvelocity_Display_1 = Show(wholegeometryvelocity_, renderView1)

#     # set active source
#     SetActiveSource(slice1)

#     # Properties modified on slice1.SliceType
#     slice1.SliceType.Origin = [X,0, 0]

#     # update the view to ensure updated data information
#     spreadSheetView1.Update()

#     # update the view to ensure updated data information
#     renderView1.Update()

#     # update the view to ensure updated data information
#     spreadSheetView2.Update()

#     # hide data in view
#     Hide(wholegeometryvelocity_, renderView1)

#     # set active source
#     SetActiveSource(calculator1)

#     # Rescale transfer function
#     calculator1Display.ScaleTransferFunction.RescaleTransferFunction(0.00444159334269, 0.307846577895)

#     # Rescale transfer function
#     calculator1Display.OpacityTransferFunction.RescaleTransferFunction(0.00444159334269, 0.307846577895)

#     # set active view
#     SetActiveView(spreadSheetView2)

#     # set active source
#     SetActiveSource(wholegeometryvelocity_)

#     # show data in view
#     cellSize1Display = Show(cellSize1, spreadSheetView2)

#     # trace defaults for the display properties.
#     cellSize1Display.FieldAssociation = 'Cell Data'

#     # set active source
#     SetActiveSource(cellSize1)

#     # save data
#     SaveData(OutputDirectory+'UpperSize.txt', proxy=cellSize1, Precision=10,
#         FieldDelimiter='\t',
#         WriteTimeSteps=1,
#         FieldAssociation='Cells')

#     # show data in view
#     calculator1Display_1 = Show(calculator1, spreadSheetView2)

#     # set active source
#     SetActiveSource(calculator1)

#     # save data
#     SaveData(OutputDirectory+'UpperVelocity.txt', proxy=calculator1, Precision=10,
#         FieldDelimiter='\t',
#         WriteTimeSteps=1,
#         FieldAssociation='Cells')

#     # set active source
#     SetActiveSource(wholegeometryvelocity_)

#     # show data in view
#     wholegeometryvelocity_Display_2 = Show(wholegeometryvelocity_, spreadSheetView2)

#     # show data in view
#     wholegeometryvelocity_Display_2 = Show(wholegeometryvelocity_, spreadSheetView2)

#     # set active view
#     SetActiveView(renderView1)

#     # show data in view
#     wholegeometryvelocity_Display_1 = Show(wholegeometryvelocity_, renderView1)

#     # set active source
#     SetActiveSource(clip1)

#     # get animation scene
#     animationScene1 = GetAnimationScene()

#     animationScene1.GoToLast()

#     # set active source
#     SetActiveSource(clip1)

#     # show data in view
#     clip1Display_1 = Show(clip1, renderView1)

#     # show color bar/color legend
#     clip1Display_1.SetScalarBarVisibility(renderView1, True)

#     # hide data in view
#     Hide(wholegeometryvelocity_, renderView1)

#     # set active source
#     SetActiveSource(wholegeometryvelocity_)

#     # show data in view
#     wholegeometryvelocity_Display_1 = Show(wholegeometryvelocity_, renderView1)

#     # set active source
#     SetActiveSource(clip1)

#     # hide data in view
#     Hide(wholegeometryvelocity_, renderView1)

#     # set active source
#     SetActiveSource(wholegeometryvelocity_)

#     # show data in view
#     wholegeometryvelocity_Display_1 = Show(wholegeometryvelocity_, renderView1)

#     Hide(wholegeometryvelocity_, renderView1)
#     # set active source
#     SetActiveSource(clip1)

#     # Properties modified on clip1.ClipType
#     clip1.ClipType.Position = [MidPoint,0, 0]
#     clip1.ClipType.Rotation = [-1.039127524240761, -0.20971833525893685, 1.506671749905828]
#     clip1.ClipType.Scale = [0.1828527359800132, 0.18210095188736036, 2.2565686673232896]

#     # update the view to ensure updated data information
#     spreadSheetView1.Update()

#     # update the view to ensure updated data information
#     renderView1.Update()

#     # update the view to ensure updated data information
#     spreadSheetView2.Update()

#     # get color legend/bar for separate_calculator1Display_VelLUT in view renderView1
#     separate_calculator1Display_VelLUTColorBar = GetScalarBar(separate_calculator1Display_VelLUT, renderView1)

#     # reset view to fit data
#     renderView1.ResetCamera()

#     # set active source
#     SetActiveSource(slice1)


#     # Properties modified on slice1.SliceType
#     slice1.SliceType.Origin = [MidPoint, 0,0]
#     slice1.SliceType.Normal = [1.0, 0.0, 0.0]

#     # update the view to ensure updated data information
#     spreadSheetView1.Update()

#     # update the view to ensure updated data information
#     renderView1.Update()

#     # update the view to ensure updated data information
#     spreadSheetView2.Update()

#     # Rescale transfer function
#     separate_calculator1Display_VelLUT.RescaleTransferFunction(-0.300152474862, 1.81186047886)

#     # Rescale transfer function
#     separate_calculator1Display_VelPWF.RescaleTransferFunction(-0.300152474862, 1.81186047886)

#     # set active source
#     SetActiveSource(clip1)

#     # set active source
#     SetActiveSource(slice1)

#     # set active source
#     SetActiveSource(calculator1)

#     # Rescale transfer function
#     calculator1Display.ScaleTransferFunction.RescaleTransferFunction(0.0139525595427, 1.7978366626)

#     # Rescale transfer function
#     calculator1Display.OpacityTransferFunction.RescaleTransferFunction(0.0139525595427, 1.7978366626)

#     # set active source
#     SetActiveSource(cellSize1)

#     # get display properties
#     cellSize1Display_1 = GetDisplayProperties(cellSize1, view=renderView1)

#     # set active view
#     SetActiveView(spreadSheetView2)

#     # set active source
#     SetActiveSource(wholegeometryvelocity_)

#     # save data
#     SaveData(OutputDirectory+'CollapsingVesselSize.txt', proxy=wholegeometryvelocity_, Precision=10,
#         FieldDelimiter='\t',
#         WriteTimeSteps=1,
#         FieldAssociation='Cells')

#     # set active view
#     SetActiveView(renderView1)

#     # get color legend/bar for velocityLUT_1 in view renderView1
#     velocityLUT_1ColorBar = GetScalarBar(velocityLUT_1, renderView1)

    
#     # set active view
#     SetActiveView(spreadSheetView2)

#     # set active view
#     SetActiveView(renderView1)

#     # show data in view
#     wholegeometryvelocity_Display_1 = Show(wholegeometryvelocity_, renderView1)

#     # set active source
#     SetActiveSource(clip1)

#     # hide data in view
#     Hide(wholegeometryvelocity_, renderView1)


#     # Properties modified on clip1.ClipType
#     clip1.ClipType.Position = [FrontPoint, 0, 0]
#     clip1.ClipType.Rotation = [-1.0391275242407632, -0.2097183352589377, 1.5066717499058302]
#     clip1.ClipType.Scale = [0.182852735980013, 0.18210095188736036, 2.2565686673232808]

#     # update the view to ensure updated data information
#     spreadSheetView1.Update()

#     # update the view to ensure updated data information
#     renderView1.Update()

#     # update the view to ensure updated data information
#     spreadSheetView2.Update()

#     # reset view to fit data
#     renderView1.ResetCamera()

#     # Properties modified on clip1.ClipType
#     clip1.ClipType.Position = [FrontPoint, 0,0]
#     clip1.ClipType.Rotation = [-1.039127524240772, -0.20971833525893951, 1.5066717499058115]
#     clip1.ClipType.Scale = [0.18285273598001084, 0.18210095188736034, 2.256568667323265]

#     # update the view to ensure updated data information
#     spreadSheetView1.Update()

#     # update the view to ensure updated data information
#     renderView1.Update()

#     # update the view to ensure updated data information
#     spreadSheetView2.Update()

#     # set active source
#     SetActiveSource(slice1)

#     # Properties modified on slice1.SliceType
#     slice1.SliceType.Origin = [FrontPoint,0,0]

#     # update the view to ensure updated data information
#     spreadSheetView1.Update()

#     # update the view to ensure updated data information
#     renderView1.Update()

#     # update the view to ensure updated data information
#     spreadSheetView2.Update()

#     # hide data in view
#     Hide(clip1, renderView1)

#     # set active source
#     SetActiveSource(calculator1)

#     # Rescale transfer function
#     calculator1Display.ScaleTransferFunction.RescaleTransferFunction(0.0123616382781, 1.71799343393)

#     # Rescale transfer function
#     calculator1Display.OpacityTransferFunction.RescaleTransferFunction(0.0123616382781, 1.71799343393)

#     # set active view
#     SetActiveView(spreadSheetView2)

#     # set active source
#     SetActiveSource(wholegeometryvelocity_)

#     # show data in view
#     cellSize1Display = Show(cellSize1, spreadSheetView2)

#     # set active source
#     SetActiveSource(cellSize1)

#     # save data
#     SaveData(OutputDirectory+'InletSize.txt', proxy=cellSize1, Precision=10,
#         FieldDelimiter='\t',
#         FieldAssociation='Cells')

#     # show data in view
#     wholegeometryvelocity_Display_2 = Show(wholegeometryvelocity_, spreadSheetView2)

#     # set active source
#     SetActiveSource(wholegeometryvelocity_)

#     # save data
#     SaveData(OutputDirectory+'InletVelocity.txt', proxy=wholegeometryvelocity_, Precision=10,
#         FieldDelimiter='\t',
#         FieldAssociation='Cells')

#     # set active view
#     SetActiveView(renderView1)

#     #### saving camera placements for all active views

#     # current camera placement for renderView1
#     renderView1.CameraPosition = [-0.0005098518773604071, 0.00024094239703430093, 0.0006290810682294644]
#     renderView1.CameraFocalPoint = [8.5539795195483e-05, -1.251353857052327e-06, -1.1271478307649734e-05]
#     renderView1.CameraViewUp = [0.34632806095905744, 0.9375473069493228, -0.03258713587310222]
#     renderView1.CameraParallelScale = 0.00013255402411569594

# #### uncomment the following to render all views
# # RenderAllViews()
# # alternatively, if you want to write images, you can use SaveScreenshot(...).


# # destroy calculator1
# Delete(calculator1)
# del calculator1

# # destroy cellSize1
# Delete(cellSize1)
# del cellSize1

# # destroy slice1
# Delete(slice1)
# del slice1


# # destroy clip1
# Delete(clip1)
# del clip1

# # destroy wholegeometryvelocity_
# Delete(wholegeometryvelocity_)
# del wholegeometryvelocity_
