import bpy
import math

gyuru = bpy.data.objects["Gyuru"].data

for vertex in gyuru.vertices:
    print("Vertex: " + str(vertex.co))
    phi = math.atan2(vertex.co[1],vertex.co[0])    
    print("phi = " + str(phi))
    if (vertex.co[2] > 0):
        vertex.co[2] +=  1.5*math.sin(-phi*7+3.14159)
    else:
        vertex.co[2] +=  -2-0.7*math.sin(phi*7)
    
gyuru.update()
