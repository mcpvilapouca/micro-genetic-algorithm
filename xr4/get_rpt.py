# -*- coding: mbcs -*-
#
# Abaqus/Viewer Release 2018 replay file
# Internal Version: 2018_03_09-14.52.59 127140
# Run by mariavp on Tue May 24 15:44:10 2022
#
filename="cube_umat"
#
import os
ppath=os.getcwd()
#
path_name=ppath+"/"+filename+".odb"
# from driverUtils import executeOnCaeGraphicsStartup
# executeOnCaeGraphicsStartup()
#: Executing "onCaeGraphicsStartup()" in the site directory ...
from abaqus import *
from abaqusConstants import *
session.Viewport(name='Viewport: 1', origin=(0.0, 0.0), width=200.666656494141,
    height=200.614044189453)
session.viewports['Viewport: 1'].makeCurrent()
session.viewports['Viewport: 1'].maximize()
from viewerModules import *
from driverUtils import executeOnCaeStartup
executeOnCaeStartup()
o1 = session.openOdb(name=path_name)
session.viewports['Viewport: 1'].setValues(displayedObject=o1)
#: Number of Assemblies:         1
#: Number of Assembly instances: 0
#: Number of Part instances:     1
#: Number of Meshes:             1
#: Number of Element Sets:       7
#: Number of Node Sets:          7
#: Number of Steps:              1
session.viewports['Viewport: 1'].odbDisplay.display.setValues(plotState=(
    CONTOURS_ON_DEF, ))
odb = session.odbs[path_name]
xyList = xyPlot.xyDataListFromField(odb=odb, outputPosition=NODAL, variable=((
    'U', NODAL, ((COMPONENT, 'U2'), )), ('SDV_PK22', INTEGRATION_POINT), ),
    nodePick=(('PART-1-1', 1, ('[#1 ]', )), ), )
xyp = session.XYPlot('XYPlot-1')
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
curveList = session.curveSet(xyData=xyList)
chart.setValues(curvesToPlot=curveList)
session.viewports['Viewport: 1'].setValues(displayedObject=xyp)
xy1 = session.xyDataObjects['_U:U2 PI: PART-1-1 N: 1']
xy2 = session.xyDataObjects['_SDV_PK22 (Avg: 75%) PI: PART-1-1 N: 1']
xy3 = combine(1+xy1, xy2)
xyp = session.xyPlots['XYPlot-1']
chartName = xyp.charts.keys()[0]
chart = xyp.charts[chartName]
c1 = session.Curve(xyData=xy3)
chart.setValues(curvesToPlot=(c1, ), )
x0 = session.xyDataObjects['_temp_2']
session.writeXYReport(fileName='report.rpt', xyData=(x0, ))
