import bpy
import math
import bmesh
import multiprocessing

angle = 12.0*math.pi/180.0

minx_threshold = 0.0
x_threshold = 0.0
removal_threshold = 1e-10

imgobj = bpy.data.objects["Par.001"]

#bpy.context.scene.objects.active = imgobj
#bpy.ops.object.mode_set(mode='EDIT')

img = imgobj.data

def dprint(s):
    # print(s)
    pass

def vnorm2(v):
    return v[0]*v[0]+v[1]*v[1]+v[2]*v[2]

def vadd(v1,v2):
    return [v1[0]+v2[0], v1[1]+v2[1], v1[2]+v2[2]]

def vmin(v1,v2):
    return [v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2]]

def vscale(v, x):
    return [v[0]*x, v[1]*x, v[2]*x]

def vcross(a,b):
    return [ a[1] * b[2] - b[1] * a[2], a[2]*b[0] - b[2]*a[0], a[0] * b[1] - b[0]*a[1] ]

def vdot(a,b):
    return a[0]*b[0] + a[1]*b[1]+a[2]*b[2]

def listkey(l):
    return str(l)

def fromlistkey(l):
    return eval(l)

zdepth = (0,0,0)
bm = bmesh.from_edit_mesh(img)
bm.verts.ensure_lookup_table()
bmverts = {}
for key, value in enumerate(bm.verts):
    bmverts[key] = value
bmloops = {}
sampleFace = None
for key, value in enumerate(bm.faces):
    sampleFace = value
    zdepth = vadd(zdepth, value.normal)
    for key2, value2 in enumerate(value.loops): 
        bmloops[value2.index] = value2
bmedges = {}
for key, value in enumerate(bm.edges):
    bmedges[key] = value

zdepth = vscale(vscale(zdepth, 1.0/math.sqrt(vdot(zdepth, zdepth))), 1.0 / math.tan(angle))
dprint("zdepth = " + str(zdepth))

# find the outer edge loops
# an edge is an inner edge if there are loops for both directions v1 -> v2 and v2-> v1
#

vertex2loop_from = {}
loop2face = {}
loops = []

def add_loop(vertex2loop, vertex_index, loop):
    if vertex_index in vertex2loop:
        loops = vertex2loop[vertex_index]
    else:
        loops = set()
        vertex2loop[vertex_index] = loops
    loops.add(loop.index)

for loopi, loop in bmloops.items():
    add_loop(vertex2loop_from, loop.vert.index, loop)

outer_edges = set()
outer_vertices = set()

no_loops = set()

vertices_in_loops = set()

branching_vertices = set()

vertices_with_no_next_loop = set()

no_facepoint_verts = set()

class VertexNormal:

    def __init__(self, vertex_co, normal_co, vertexObject):
        self.vertex = vertex_co
        self.normal = normal_co
        self.vertexObject = vertexObject

    def __getstate__(self):
        return (self.vertex[0], self.vertex[1], self.vertex[2], self.normal[0], self.normal[1], self.normal[2])

    def __setstate__(self, d):
        self.vertex = (d[0], d[1], d[2])
        self.normal = (d[3], d[4], d[5])
        self.vertexObject = None

    def __str__(self):
        return "(vertex: " + str(self.vertex) + ", normal: " + str(self.normal) + ")"

    def __repr__(self):
        return "(vertex: " + str(self.vertex) + ", normal: " + str(self.normal) + ")"

def get_normal(pv, v, nv):
    nvv = vmin(nv,v)
    vpv = vmin(v,pv)
    n1 = vcross(zdepth, nvv)
    n2 = vcross(zdepth, vpv)
    if vnorm2(n1) > 1e-4 * vnorm2(zdepth) * vnorm2(nvv):
        n1 = vscale(n1, 1.0/math.sqrt(vnorm2(n1)))
    else:
        n1 = (0.,0.,0.)
    if vnorm2(n2) > 1e-4 * vnorm2(zdepth) * vnorm2(vpv):
        n2 = vscale(n2, 1.0/math.sqrt(vnorm2(n2)))
    else:
        n2 = (0.,0.,0.)
    n = vadd(n1, n2)
    if vnorm2(n) > 1e-10:
        n = vscale(n, 1.0/math.sqrt(vnorm2(n)))
        if 1. + vdot(n1, n2) > 1e-2:
            n = vscale(n, 1.0/ math.sqrt(0.5*(1.+vdot(n1,n2))))
        else:
            n = (0., 0., 0.)
    else:
        n = (0.,0.,0.)
    return n

def test_removal(polyindex1, poly1, i1):
    # next index
    ni1 = i1+1
    if ni1 == len(poly1):
        ni1 = 0
    pi1 = i1-1
    if pi1 < 0:
        pi1 = len(poly1)-1
    v1 = poly1[i1].vertex
    nv1 = poly1[ni1].vertex
    n1 = poly1[i1].normal
    nn1 = poly1[ni1].normal
    pv1 = poly1[pi1].vertex
    pn1 = poly1[pi1].normal
    #
    # calculate removals
    #
    # v1 + x n1 - nv1 - x nn1 = u * (pv1 + x pn1 - v1 - x n1)
    #
    # dotprod(z, crossprod(v1 + x n1 - nv1 - x nn1, pv1 + x pn1 - v1 - x n1)) = 0
    # -> 2nd order equation for x.
    #
    # a(b x c) = a_i b_j c_k epsilon_ijk = (a x b) c
    #
    # dotprod( crossprod(z , v1 - nv1) + x* crossprod(z, n1 - nn1), pv1  - v1 + x * (pn1 - n1)) = 0
    #
    # x^2 * dotprod( crossprod(z, n1 - nn1), pn1 - n1)
    # + x* ( dotprod( crossprod(z , v1 - nv1), pn1  - n1)  + dotprod(crossprod(z, n1 - nn1), pv1 - v1))
    # + dotprod( crossprod(z , v1 - nv1), pv1  - v1)  = 0
    #
    # for all solutions calculate sign of 'u' by sgn(dotprod(v1 + x n1 - nv1 - x nn1, pv1 + x pn1 - v1 - x n1)), if it is -1 then x is a good solution.
    # 
    #
    a = vdot(vcross(zdepth, vmin(n1, nn1)), vmin(pn1, n1))
    b = vdot(vcross(zdepth, vmin(v1, nv1)), vmin(pn1, n1)) + vdot(vcross(zdepth, vmin(n1, nn1)), vmin(pv1, v1))
    c = vdot(vcross(zdepth, vmin(v1, nv1)), vmin(pv1, v1))
    if math.fabs(a) <= 1e-10*(math.fabs(c) + math.fabs(b)):
        if math.fabs(b) <= 1e-10*math.fabs(c):
            x1 = None
            x2 = None
        else:
            x1 = -c / b
            x2 = None
    elif b*b - 4*a*c >= 0:
        x1 = (-b + math.sqrt(b*b - 4*a*c))/(2*a)
        x2 = (-b - math.sqrt(b*b - 4*a*c))/(2*a)
    else:
        x1 = None
        x2 = None
    dprint("removal test " + str([polyindex1, i1]) + ": x1 = " + str(x1) + ", x2 = " + str(x2))
    if not(x1 is None):
        if x1 > 0:
            va = vadd(vmin(v1, nv1), vscale(vmin(n1, nn1), x1))
            vb = vadd(vmin(pv1, v1), vscale(vmin(pn1, n1), x1))
            if vnorm2(va) < removal_threshold * vnorm2(vb) or vnorm2(vb) < removal_threshold * vnorm2(va):
                dprint("removal test " + str([polyindex1, i1]) + ": vnorm2(va) = " + str(vnorm2(va)) + ", vnorm2(vb) = " + str(vnorm2(vb)))
                pass
            else:
                u = vdot(va, vb)
                dprint("removal test " + str([polyindex1, i1]) + ": u = " + str(u) + ", va = " + str(va) + ", vb = " + str(vb))
                if u >= 0:
                    x1 = None
        else:
            x1 = None
    if not(x2 is None):
        if x2 > 0:
            va = vadd(vmin(v1, nv1), vscale(vmin(n1, nn1), x2))
            vb = vadd(vmin(pv1, v1), vscale(vmin(pn1, n1), x2))
            if vnorm2(va) < removal_threshold * vnorm2(vb) or vnorm2(vb) < removal_threshold * vnorm2(va):
                dprint("removal test " + str([polyindex1, i1]) + ": vnorm2(va) = " + str(vnorm2(va)) + ", vnorm2(vb) = " + str(vnorm2(vb)))
                pass
            else:
                u = vdot(va, vb)
                dprint("removal test " + str([polyindex1, i1]) + ": u = " + str(u) + ", va = " + str(va) + ", vb = " + str(vb))
                if u >= 0:
                    x2 = None
        else:
            x2 = None
    if not(x1 is None) and not(x2 is None):
        if x1 < x2:
            x2 = None
        else:
            x1 = x2
            x2 = None
    x_result = None
    if not(x1 is None):
        x_result = x1
    return x_result


def test_collision(polyindex1, poly1, i1, polyindex2, poly2, i2):
    result = [None, None]
    ni1 = i1+1
    if ni1 == len(poly1):
        ni1 = 0
    v1 = poly1[i1].vertex
    nv1 = poly1[ni1].vertex
    n1 = poly1[i1].normal
    nn1 = poly1[ni1].normal
    v2 = poly2[i2].vertex
    n2 = poly2[i2].normal
    #
    # v1 + x*n1  + alpha * (nv1 + x*nn1 - v1 - x*n1) = v2 + x*n2
    # 0 <= alpha <= 1
    #
    #
    # v1-v2 + x*(n1-n2) = alpha* (v1-nv1 + x*(n1-nn1))
    #
    # if v1 - v2 || v1-nv1 --> we are ready, no overlap
    # else dot product by these:
    #
    # dot(v1-v2, v1-v2) + x*dot(n1-n2,v1-v2) = alpha*(dot(v1-nv1, v1-v2) + x*dot(n1-nn1, v1-v2))
    # dot(v1-v2, v1-nv1) + x*dot(n1-n2,v1-nv1) = alpha*(dot(v1-nv1, v1-nv1) + x*dot(n1-nn1, v1-nv1))
    #
    # (dot(v1-v2, v1-nv1) + x*dot(n1-n2,v1-nv1))*(dot(v1-nv1, v1-v2) + x*dot(n1-nn1, v1-v2)) = (dot(v1-v2, v1-v2) + x*dot(n1-n2,v1-v2))*(dot(v1-nv1, v1-nv1) + x*dot(n1-nn1, v1-nv1))
    #
    # x^2* ( dot(n1-n2,v1-nv1)* dot(n1-nn1, v1-v2) - dot(n1-n2,v1-v2)*dot(n1-nn1, v1-nv1) ) + x *( dot(n1-n2,v1-nv1) * dot(v1-nv1, v1-v2) + dot(v1-v2, v1-nv1) * dot(n1-nn1,v1-v2) - dot(n1-n2,v1-v2) * dot(v1-nv1, v1-nv1) - dot(v1-v2, v1-v2) * dot(n1-nn1, v1-nv1)) + dot(v1-v2, v1-nv1) * dot(v1-nv1, v1-v2) - dot(v1-v2, v1-v2) * dot(v1-nv1, v1-nv1) = 0
    #
    #
    a = vdot(vmin(n1,n2),vmin(v1,nv1))*vdot(vmin(n1,nn1), vmin(v1,v2)) - vdot(vmin(n1,n2),vmin(v1,v2))*vdot(vmin(n1,nn1),vmin(v1,nv1))
    b = vdot(vmin(n1,n2),vmin(v1,nv1))*vdot(vmin(v1,nv1), vmin(v1,v2)) + vdot(vmin(v1,v2),vmin(v1,nv1))*vdot(vmin(n1,nn1),vmin(v1,v2)) - vdot(vmin(n1,n2),vmin(v1,v2))*vdot(vmin(v1,nv1),vmin(v1,nv1))-vdot(vmin(v1,v2),vmin(v1,v2))*vdot(vmin(n1,nn1),vmin(v1,nv1))
    c = vdot(vmin(v1,v2), vmin(v1,nv1))*vdot(vmin(v1,nv1), vmin(v1,v2)) - vdot(vmin(v1,v2), vmin(v1,v2))*vdot(vmin(v1,nv1), vmin(v1,nv1))
    d = b*b - 4 * a *c
    if a == 0.0:
        if b == 0.0:
            x1 = -1
            x2 = -1
        else:
            x1 = -c/b
            x2 = -c/b
    elif d >= 0.0:
        x1 = (-b + math.sqrt(d))/2.0/a
        x2 = (-b - math.sqrt(d))/2.0/a
    else:
        x1 = -1
        x2 = -1
        # verify result
    divisor1 = vdot(vmin(v1,nv1), vmin(v1,nv1)) + x1*vdot(vmin(n1,nn1), vmin(v1,nv1))
    if x1 > x_threshold and math.fabs(divisor1) > 0.0:
        alpha1 = (vdot(vmin(v1,v2), vmin(v1,nv1)) + x1*vdot(vmin(n1,n2),vmin(v1,nv1))) /  divisor1
        dprint("addition test 1: " + str([polyindex1, i1, ni1, polyindex2, i2, x1, alpha1]))
        # v1 + x*n1  + alpha * (nv1 + x*nn1 - v1 - x*n1) = v2 + x*n2
        v1pxn1 = vadd(vscale(n1, x1), v1)
        nv1pxnn1 = vadd(vscale(nn1, x1), nv1)
        left = vadd(vscale(v1pxn1,1.0-alpha1), vscale(nv1pxnn1, alpha1))
        right = vadd(v2, vscale(n2, x1))
        diff = vmin(left, right)
        if vnorm2(diff) > 1e-10:
            dprint("Equation not solved correctly (1). v1 = " + str(v1) + ", n1 = " + str(n1) + ", nv1 = " + str(nv1) + ", nn1 = " + str(nn1) + ", v2 = " + str(v2) + ", n2 = " + str(n2) + ", left = " + str(left) + ", right = " + str(right) + ", x = " + str(x1) + ", alpha = " + str(alpha1) + ", v1+x n1 = " + str(v1pxn1) + ", nv1 + x nn1 = " + str(nv1pxnn1))
        if alpha1 >= 0.0 and alpha1 <= 1.0:
            result = [x1, alpha1]
            if vnorm2(diff) > 1e-10:
                print("Equation not solved correctly (1). v1 = " + str(v1) + ", n1 = " + str(n1) + ", nv1 = " + str(nv1) + ", nn1 = " + str(nn1) + ", v2 = " + str(v2) + ", n2 = " + str(n2) + ", left = " + str(left) + ", right = " + str(right) + ", x = " + str(x1) + ", alpha = " + str(alpha1) + ", v1+x n1 = " + str(v1pxn1) + ", nv1 + x nn1 = " + str(nv1pxnn1))
    divisor2 = vdot(vmin(v1,nv1), vmin(v1,nv1)) + x2*vdot(vmin(n1,nn1), vmin(v1,nv1))
    if x2 > x_threshold and  math.fabs(divisor2) > 0.0:

        alpha2 = (vdot(vmin(v1,v2), vmin(v1,nv1)) + x2*vdot(vmin(n1,n2),vmin(v1,nv1))) /  (vdot(vmin(v1,nv1), vmin(v1,nv1)) + x2*vdot(vmin(n1,nn1), vmin(v1,nv1)))
        dprint("addition test 2: " + str([polyindex1, i1, ni1, polyindex2, i2, x2, alpha2]))
        # v1 + x*n1  + alpha * (nv1 + x*nn1 - v1 - x*n1) = v2 + x*n2
        v1pxn1 = vadd(vscale(n1, x2), v1)
        nv1pxnn1 = vadd(vscale(nn1, x2), nv1)
        left = vadd(vscale(v1pxn1,1.0-alpha2), vscale(nv1pxnn1, alpha2))
        right = vadd(v2, vscale(n2, x2))
        diff = vmin(left, right)
        if vnorm2(diff) > 1e-10:
            dprint("Equation not solved correctly (2). v1 = " + str(v1) + ", n1 = " + str(n1) + ", nv1 = " + str(nv1) + ", nn1 = " + str(nn1) + ", v2 = " + str(v2) + ", n2 = " + str(n2) + ", left = " + str(left) + ", right = " + str(right) + ", x = " + str(x2) + ", alpha = " + str(alpha2) + ", v1+x n1 = " + str(v1pxn1) + ", nv1 + x nn1 = " + str(nv1pxnn1))
        if alpha2 >= 0.0 and alpha2 <= 1.0 and vnorm2(diff) > 1e-10:
            print("Equation not solved correctly (2). v1 = " + str(v1) + ", n1 = " + str(n1) + ", nv1 = " + str(nv1) + ", nn1 = " + str(nn1) + ", v2 = " + str(v2) + ", n2 = " + str(n2) + ", left = " + str(left) + ", right = " + str(right) + ", x = " + str(x2) + ", alpha = " + str(alpha2) + ", v1+x n1 = " + str(v1pxn1) + ", nv1 + x nn1 = " + str(nv1pxnn1))
    if x2 > x_threshold and (x1 < x_threshold or x2 < x1) and alpha2 >= 0.0 and alpha2 <= 1.0:
        result = [x2, alpha2]
    return result

def call_test_removal(package):
    results = []
    for w in package:
        results.append(test_removal(w[0], w[1], w[2]))
    return results

def get_test_removal(results, package, removal_result):
    for i in range(0, len(package)):
        w = package[i]
        removal_result[listkey([w[0], w[2]])] = results[i]

def call_test_collision(package):
    results = []
    for w in package:
        results.append(test_collision(w[0], w[1], w[2], w[3], w[4], w[5]))
    return results

def get_test_collision(results, package, collision_result):
    for i in range(0, len(package)):
        w = package[i]
        collision_result[listkey([w[0], w[2], w[3], w[5]])] = results[i]


def max_shift(loopvertices, collisions, removals):
    # two types of collisions exist:
    # 1) a shifted vertex collides with a shifted edge
    #    v' + alpha*(nv'-v') = v2'
    #    remediation: introduce a middle point in v' - nv' with zero normal. lock v2 to that point with zero normal as well.
    # 2) adjacent vertices cross
    #    pv' - v' - nv' in one line, and pv' - v' is opposite to v' - nv'.
    #
    #    v' - nv' = u * (pv' - v')
    #    u < 0
    #    remediation:
    #    remove v'
    #
    collisions.clear()
    removals.clear()
    removal_result = {}
    collision_result = {}
    minx = 10000.0
    cpus = multiprocessing.cpu_count()
    packages=[]
    for t in range(0, cpus):
        packages.append([])
    cpuindex=0
    for polyindex1 in range(0, len(loopvertices)):
        poly1 = loopvertices[polyindex1]
        for i1 in range(0, len(poly1)):
            package = packages[cpuindex]
            cpuindex = (cpuindex + 1) % cpus
            package.append((polyindex1, poly1, i1))
    results = []
    pool = multiprocessing.Pool(cpus)
    for t in range(0, cpus):
        results.append(pool.apply_async(call_test_removal, [packages[t]]))
    pool.close()
    pool.join()
    for t in range(0, cpus):
        get_test_removal(results[t].get(), packages[t], removal_result)
        packages[t] = []
    for polyindex1 in range(0, len(loopvertices)):
        poly1 = loopvertices[polyindex1]
        for i1 in range(0, len(poly1)):
            #
            # test if vertex collides with another edge.
            #
            ni1 = i1+1
            if ni1 == len(poly1):
                ni1 = 0
            for polyindex2 in range(0, len(loopvertices)):
                poly2 = loopvertices[polyindex2]
                for i2 in range(0, len(poly2)):
                    if (polyindex1 != polyindex2 or (i1 != i2 and ni1 != i2)):
                        package = packages[cpuindex]
                        cpuindex = (cpuindex + 1) % cpus
                        package.append((polyindex1, poly1, i1, polyindex2, poly2, i2))
    results = []
    pool = multiprocessing.Pool(cpus)
    for t in range(0, cpus):
        results.append(pool.apply_async(call_test_collision, [packages[t]]))
    pool.close()
    pool.join()
    for t in range(0, cpus):
        get_test_collision(results[t].get(), packages[t], collision_result)
    for polyindex1 in range(0, len(loopvertices)):
        poly1 = loopvertices[polyindex1]
        for i1 in range(0, len(poly1)):
            x = removal_result[listkey([polyindex1, i1])]
            if not(x is None):
                if x <= minx + minx_threshold:
                    if x < minx - minx_threshold:
                        collisions.clear()
                        removals.clear()
                        minx = x
                    removals[listkey([polyindex1, i1])] = [1, 0]
    for polyindex1 in range(0, len(loopvertices)):
        poly1 = loopvertices[polyindex1]
        for i1 in range(0, len(poly1)):
            #
            # test if vertex collides with another edge.
            #
            ni1 = i1+1
            if ni1 == len(poly1):
                ni1 = 0
            for polyindex2 in range(0, len(loopvertices)):
                poly2 = loopvertices[polyindex2]
                for i2 in range(0, len(poly2)):
                    if (polyindex1 != polyindex2 or (i1 != i2 and ni1 != i2)) and not(listkey([polyindex2, i2]) in removals):
                        [x, alpha] = collision_result[listkey([polyindex1, i1, polyindex2, i2])]
                        if not(x is None) and x <= minx + minx_threshold:
                            if x < minx - minx_threshold:
                                collisions.clear()
                                removals.clear()
                                minx = x
                            if not(listkey([polyindex1, i1]) in removals) and not(listkey([polyindex1, ni1]) in removals):
                                collisions[listkey([polyindex2, i2])] = [polyindex1, i1, ni1, alpha]
    if minx == 10000.0:
        minx = -1
    else:
        if len(removals) > 0 and len(collisions) > 0:
            dprint("Both removals and collisions!, removals = " + str(removals) + ", collisions = " + str(collisions))
        if len(removals) > 1:
            dprint("Multi-removals: " + str(removals))
        if len(collisions) > 1:
            dprint("Multi-collisions: " + str(collisions))
    print("minx = " + str(minx) + ", collisions = " + str(collisions) + ", removals = " + str(removals))
    return minx

def addface(bm, verts, sampleFace):
    dprint("addface: verts = " + str(verts[0].co) + ", " + str(verts[1].co) + ", " + str(verts[2].co))
    # face normal points to outside:
    newface1 = bm.faces.new( verts, sampleFace)
    newface1.normal_update()
    if vdot(newface1.normal, zdepth) > math.sin(angle)*0.5*math.sqrt(vnorm2(zdepth)*vnorm2(newface1.normal)):
        newface1.normal_flip()
    elif vdot(newface1.normal, zdepth) > -math.sin(angle)*0.5*math.sqrt(vnorm2(zdepth)*vnorm2(newface1.normal)):
        dprint("Guessing right normal: " + str(newface1.normal) + ", face: " + str(newface1))
        direction = (0.,0.,0.)
        maxdot = 0
        for ek,e in enumerate(newface1.edges):
            for otherfacekey,otherface in enumerate(e.link_faces):
                if not(otherface is newface1):
                    facedot = math.fabs(vdot(otherface.normal, newface1.normal))
                    if facedot > maxdot:
                        maxdot = facedot
                        direction = otherface.normal
                        dprint("otherface: " + str(otherface))
        if maxdot > 0:
            dprint("maxdot = " + str(maxdot) + ", direction = " + str(direction))
            if vdot(newface1.normal, direction) < 0:
                dprint("flip!")
                newface1.normal_flip()
            else:
                dprint("no flip!")
        else:
            dprint("Couldn't decide direction of the new face")
    dprint("added face: " + str(newface1) + ", normal: " + str(newface1.normal))
    
for edgei, edge in bmedges.items():
    from_found = False
    to_found = False
    if not (edge.verts[0].index in vertex2loop_from):
        print("Vertex doesn't have a loop starting from it: " + str(edge.verts[0].index))
        no_loops.add(edge.verts[0].index)
    else:
        for loop in vertex2loop_from[edge.verts[0].index]:
            loop_edge = bmloops.get(loop).edge
            if loop_edge.verts[0].index == edge.verts[1].index or loop_edge.verts[1].index == edge.verts[1].index:
                from_found = True
    if not (edge.verts[1].index in vertex2loop_from):
        print("Vertex doesn't have a loop starting from it: " + str(edge.verts[1].index))
        no_loops.add(edge.verts[1].index)
    else:
        for loop in vertex2loop_from[edge.verts[1].index]:
            loop_edge = bmloops[loop].edge
            if loop_edge.verts[0].index == edge.verts[0].index or loop_edge.verts[1].index == edge.verts[0].index:
                to_found = True
    if (to_found or from_found) and not (to_found and from_found):
        outer_edges.add(edge.index)
        outer_vertices.add(edge.verts[0].index)
        outer_vertices.add(edge.verts[1].index)

if len(outer_vertices) != len(bmverts):
    print("There are inner vertices, outer_vertices size = " + str(len(outer_vertices)) + ", all vertices count = " + str(len(bmverts)))
    bm.select_flush(True)
    bm.select_mode = {"VERT"}
    for vertexi, vertex in bmverts.items():
        if not (vertex.index in outer_vertices):
            print("Inner vertex: " + str(vertex.index) + ": " + str(vertex.co[0]) + ", " + str(vertex.co[1]) + ", " + str(vertex.co[2]))
            vertex.select = True
        else:
            vertex.select = False
elif len(no_loops) > 0:
    bm.select_flush(True)
    bm.select_mode = {"VERT"}
    for vertexi, vertex in bmverts.items():
        vertex.select = (vertex.index in no_loops)
    print("Vertices without loops from them selected.")
else:
    print("Good. There are no inner vertices.")
        

    
    for loop_vertex in outer_vertices:
        if not (loop_vertex in vertices_in_loops):
            dprint("Starting from vertex " + str(loop_vertex))
            loop = []
            closed = False
            while not closed:
                next_loops = set()
                if not loop_vertex in vertex2loop_from:
                    dprint("No loop from vertex " + str(loop_vertex))
                    break
                for loop_candidate in vertex2loop_from[loop_vertex]:
                    dprint("loop candidate: " + str(loop_candidate))
                    if bmloops[loop_candidate].edge.index in outer_edges:
                        dprint("in outer edges")
                        if bmloops[loop_candidate].edge.verts[0].index == loop_vertex:
                            to_vertex = bmloops[loop_candidate].edge.verts[1].index
                        else:
                            to_vertex = bmloops[loop_candidate].edge.verts[0].index
                        dprint("to_vertex = " + str(to_vertex))
                        if len(loop) > 0:
                            first_vertex = bmloops[loop[0]].vert.index
                        else:
                            first_vertex = -1
                        if not (to_vertex in vertices_in_loops) or to_vertex == first_vertex:
                            duplicate = False
                            for l in next_loops:
                                if bmloops[loop_candidate].edge.index == bmloops[l].edge.index:
                                    print("Duplicate loop: " + str(loop_candidate) + " and " + str(l))
                                    duplicate = True
                            if not duplicate:
                                dprint("adding to next_loops")
                                next_loops.add(loop_candidate)
                if len(next_loops) != 1:
                    if len(next_loops) > 1:
                        branching_vertices.add(loop_vertex)
                    for l in next_loops:
                        dprint("loop: " + str(l) + ": vertex_index = " + str(bmloops[l].vert.index) + ", edge: " + str(bmloops[l].edge.index) + ", edge[0] = " + str(bmloops[l].edge.verts[0].index) + ", edge[1] = " + str(bmloops[l].edge.verts[1].index))
                    dprint("Next loop couldn't be found for vertex " + str(loop_vertex) + " found: " + str(len(next_loops)) + ": " + str(next_loops))
                    if len(loop) > 1:
                        dprint("loop: " + str(loop))
                        vertices_with_no_next_loop.add(loop_vertex)
                    break
                else:
                    for loop_candidate in next_loops:
                        dprint("Found next loop: " + str(loop_candidate) + ", from vertex: " + str(loop_vertex))
                        vertices_in_loops.add(loop_vertex)
                        dprint(" to vertex " + str(to_vertex))
                        loop_vertex = to_vertex
                        loop.append(loop_candidate)
                        dprint("First vertex: " + str(bmloops[loop[0]].vert.index) + ", edge.vertices[0] = " + str(bmloops[loop[0]].edge.verts[0].index) + ", edge.vertices[1] = " + str(bmloops[loop[0]].edge.verts[1].index))
                        if to_vertex == bmloops[loop[0]].vert.index:
                            closed = True
                            loops.append(loop)
                            dprint("closed")
    dprint("Found " + str(len(loops)) + " loops.")

if len(branching_vertices) > 0:
    bm.select_flush(True)
    bm.select_mode = {"VERT"}
    for vertex in branching_vertices:
        bmverts[vertex].select = True
    print("Problem: branching vertices in outer edge loops found and selected")
elif len(vertices_with_no_next_loop) > 0:
    bm.select_flush(True)
    bm.select_mode = {"VERT"}
    for vertex in vertices_with_no_next_loop:
        bmverts[vertex].select = True
    print("Problem: vertices with no next loop in outer edge loops found and selected")
else:
    #  check if all vertices are selected
    if len(vertices_in_loops) != len(outer_vertices):
        print("Not all outer vertices have been found while checking the loops")
    else:
        print("So far so good, selecting outer edge.")
        bm.select_flush(True)
        bm.select_mode = {"EDGE"}
        for edge in outer_edges:
            bmedges[edge].select = True
        
        # generate the normal vector for each vertex
        polynormals = []
        for polyindex in range(0, len(loops)):
            poly = loops[polyindex]
            normals = []
            for i in range(0, len(poly)):
                # previous index
                pi = i-1
                if pi < 0:
                    pi = len(poly)-1
                # next index
                ni = i+1
                if ni == len(poly):
                    ni = 0
                pv = bmloops[poly[pi]].vert.co
                v = bmloops[poly[i]].vert.co
                nv = bmloops[poly[ni]].vert.co
                n = get_normal(pv, v, nv)
                if math.fabs(vdot(n, zdepth)) > 1e-8:
                    dprint("Incorrect normal: " + str(i) + ", v = " + str(v) + ", pv = " + str(pv) + ", nv = " + str(nv) + ", zdepth = " + str(zdepth) + ", nv - v = " + str(nvv) + ", v - pv = " + str(vpv) + ", n1 = " +str(n1) + ", n2 = " + str(n2) + ", n = " + str(n) + ", z x (nv -v) = " + str(vcross(zdepth, nvv)) + ", z x (v - pv) = " + str(vcross(zdepth, vpv)))
                #                dprint("n = " + str(n) + ", nvv = " + str(nvv) + ", vpv = " + str(vpv))
                #sum_faces = 0.0
                #sum_faces_n = 0.0
                #count_faces = 0.0
                #for fi, f in enumerate(bmloops[poly[i]].vert.link_faces):
                #    sum_faces += f.calc_area()
                #    count_faces += 1.0
                #v_shift = vadd(v, vscale(n, math.sqrt(sum_faces/count_faces)*0.01))
                #for fi, f in enumerate(bmloops[poly[i]].vert.link_faces):
                #    shifted = False
                #    shift_index = -1
                #    for fvi, fv in enumerate(f.verts):
                #        if vnorm2(vmin(fv.co,v)) == 0.0:
                #            fv.co = v_shift
                #            shift_index = fvi
                #            shifted = True
                #    if not shifted:
                #        dprint("ERROR: vector to shift not found in copied face.")
                #sum_faces_n += f.calc_area()
                # shift back
                #    for fvi, fv in enumerate(f.verts):
                #        if fvi == shift_index:
                #            fv.co = v
                #dprint("sum_faces = " + str(sum_faces) + ", sum_faces_n = " + str(sum_faces_n))
                #if sum_faces_n > sum_faces:
                #    n = vscale(n, -1.)
                #    dprint("norm flip needed")
                normals.append(n)
                # dprint("normal[" + str(i) + "] = " + str(n) + " norm2 = " + str(vnorm2(n)) + " nv-v_n = " + str(nvvn) + " v-pv_n = " + str(vpvn) + ", x = " + str(x))
            polynormals.append(normals)
            dprint("A polygon has been processed: " + str(len(normals)) + " normals, loops: " + str(len(poly)))
        dprint("Normals have been calculated: " + str(len(polynormals)))
        if len(no_facepoint_verts) > 0:
            bm.select_flush(True)
            bm.select_mode = {"VERT"}
            for vi, v in bmverts.items():
                v.select = False
            for v in no_facepoint_verts:
                bmverts[v].select = True
            dprint("Vertices with no face point outside of the edge found and selected.")
        else:
            # data is transformed to keep only vertices
            loopvertices = []
            # we also keep track of faces we have to add for closing the 'top' - doesn't work well due to intersecting faces
            # closingfaces = {}
            # vertindex = {}
            for polyindex in range(0, len(loops)):
                vertices = []
                loopvertices.append(vertices)
                poly = loops[polyindex]
                for i in range(0, len(poly)):
                    normal_co = polynormals[polyindex][i]
                    vertex_co = bmloops[poly[i]].vert.co
                    # vertindex[bmloops[poly[i]].vert.index] = [polyindex, i]
                    vertices.append(VertexNormal(vertex_co, normal_co, bmloops[poly[i]].vert))
            # for polyindex in range(0, len(loops)):
            #     poly = loops[polyindex]
            #     for i in range(0, len(poly)):
            #         for fk, f in enumerate(bmloops[poly[i]].vert.link_faces):
            #             face = []
            #             for vk, v in enumerate(f.verts):
            #                 if v.index in vertindex:
            #                     face.append(vertindex[v.index])
            #                 else:
            #                     print("Vertex " + str(v) + " not found in index, v.index = " + str(v.index))
            #             if not (str(face) in closingfaces):
            #                 closingfaces[str(face)] = face
            # dprint("Closing faces = " + str(closingfaces))
            # how much we can shift without overlapping is calculated
            # we shift until no shift is possible
            # collisions will be selected
            bm.select_flush(True)
            bm.select_mode = {"VERT"}
            for vertexi, vertex in bmverts.items():
                vertex.select = False
            iter = 1
            while True:
                collisions = {}
                removals = {}
                # recalculate all normals

                for polyindex in range(0, len(loopvertices)):
                    poly = loopvertices[polyindex]
                    for i in range(0, len(poly)):
                        ni = i+1
                        if ni == len(poly):
                            ni = 0
                        pi = i-1
                        if pi < 0:
                            pi = len(poly)-1
                    
                        poly[i].normal = get_normal(poly[pi].vertex, poly[i].vertex, poly[ni].vertex)

                x = max_shift(loopvertices, collisions, removals)
                if x < 0:
                    print("No vertices can be shifted, stopping")
                    break
                dprint("Vertex shift amount: " + str(x) + ", collisions: " + str(collisions) + ", removals: " + str(removals))
                # create a set of new loops from the old loops:
                # shift all vertices
                newloopvertices = []
                additions = {}
                loopfaces = []
                for polyindex in range(0, len(loopvertices)):
                    newpoly = []
                    faces = []
                    loopfaces.append(faces)
                    newloopvertices.append(newpoly)
                    poly = loopvertices[polyindex]
                    for i in range(0, len(poly)):
                        ni = i+1
                        if ni == len(poly):
                            ni = 0
                        pi = i-1
                        if pi < 0:
                            pi = len(poly)-1
                        newvert = bm.verts.new(vadd(vadd(vscale(poly[i].normal, x),poly[i].vertex), vscale(zdepth, x)))
                        newVertexNormal = VertexNormal(newvert.co, poly[i].normal, newvert)
                        dprint("Loop " + str(polyindex) + ", new vertex at " + str(len(newpoly)) + ": " + str(newvert.co))
                        newpoly.append(newVertexNormal)
                        faces.append([0, i, 0, ni, 1, i])
                        dprint("Loop " + str(polyindex) + ", Added face: " + str(faces[-1]))
                        faces.append([0, ni, 1, ni, 1, i])
                        dprint("Loop " + str(polyindex) + ", Added face: " + str(faces[-1]))
                        if listkey([polyindex, i]) in collisions:
                            #  collisions[listkey([polyindex2, i2])] = [polyindex1, i1, ni1, alpha]
                            collision = collisions[listkey([polyindex, i])]
                            cpi = collision[1]-1
                            if cpi < 0:
                                cpi = len(poly)-1
                            cnni = collision[2] + 1
                            if cnni == len(poly):
                                cnni = 0
                            add = True
                            if collision[0] == polyindex and i == cpi:
                                removals[listkey([polyindex, collision[1]])] = [3,0]
                                dprint("Loop " + str(polyindex) + ", Removal: " + str([polyindex, collision[1]]))
                                add = False
                            if collision[0] == polyindex and i == cnni:
                                removals[listkey([polyindex, collision[2]])] = [4,0]
                                dprint("Loop " + str(polyindex) + ", Removal: " + str([polyindex, collision[2]]))
                                add = False
                            newvert.select = True
                            if add:
                                # two edges intersect, and they don't share a common vertex
                                #
                                additions[listkey([collision[0], collision[1]])] = [polyindex, i, collision[3]]
                        else:
                            if listkey([polyindex, i]) in removals:
                                dprint("Loop " + str(polyindex) + ", Removal: " + str(i))
                                poly[i].vertexObject.select = True
                            else:
                                newvert.select = False
                # process additions
                #
                # If the loop intersects with itself: generate two new loops.
                # If two loops intersect: join the loops into a single one.
                #
                orig_additions = {}
                for k,v in additions.items():
                    orig_additions[k] = v
                orig_removals = {}
                for k,v in removals.items():
                    orig_removals[k] = v
                if len(additions) > 0:
                    # add faces first. We will start anew. Removals are not processed.
                    if len(removals) > 0:
                        print("ERROR: removals found but additions as well. Not processing removals.");
                        
                    for polyindex in range(0, len(loopvertices)): 
                        faces = loopfaces[polyindex]
                        poly = loopvertices[polyindex]
                        newpoly = newloopvertices[polyindex]
                        for face in faces:
                            dprint("Loop " + str(polyindex) + ", Face creation: " + str(face))                            
                            if face[2] == 0:
                                addface( bm, (poly[face[1]].vertexObject, poly[face[3]].vertexObject, newpoly[face[5]].vertexObject), sampleFace)
                            else:
                                if face[3] != face[5]:
                                    if listkey([polyindex,face[5]]) in additions:
                                        # add a split face if we will insert a new vertex when continuing
                                        addition = additions[listkey([polyindex,face[5]])]
                                        # the collision point:
                                        midVertex = newloopvertices[addition[0]][addition[1]].vertexObject
                                        addface( bm, (poly[face[1]].vertexObject, newpoly[face[3]].vertexObject, midVertex), sampleFace)
                                        addface( bm, (poly[face[1]].vertexObject, midVertex, newpoly[face[5]].vertexObject), sampleFace)
                                    else:
                                        addface( bm, (poly[face[1]].vertexObject, newpoly[face[3]].vertexObject, newpoly[face[5]].vertexObject), sampleFace)
                    while len(additions) > 0:
                        for keya, va in additions.items():
                            a = fromlistkey(keya)
                            if a[0] == va[0]:
                                # self-collision, generate two loops
                                # collision point is a[1]*(1 - va[2]) + nextof(a[1])*(va[2]) = va[1]
                                #
                                # second loop will be nextof(a[1])...va[1], first loop will be va[1]...a[1].
                                #
                                na1 = a[1] + 1
                                if na1 == len(newloopvertices[a[0]]):
                                    na1 = 0
                                newpoly1=[]
                                newpoly2=[]
                                newpoly2_indexes = {}
                                newpoly1_indexes = {}
                                i = na1
                                while i != va[1]:
                                    newpoly1_indexes[i] = len(newpoly1)
                                    newpoly1.append(newloopvertices[a[0]][i])
                                    i = i + 1
                                    if i == len(newloopvertices[a[0]]):
                                        i = 0
                                newpoly1_indexes[i] = len(newpoly1)
                                newVertexNormal = VertexNormal(newloopvertices[a[0]][i].vertex, (0.0, 0.0, 0.0), newloopvertices[a[0]][i].vertexObject)
                                newpoly1.append(newVertexNormal)
                                while i != na1:
                                    newpoly2_indexes[i] = len(newpoly2)
                                    newpoly2.append(newloopvertices[a[0]][i])
                                    i = i + 1
                                    if i == len(newloopvertices[a[0]]):
                                        i = 0
                                newpoly2[0].normal = get_normal(newpoly2[len(newpoly2)-1].vertex, newpoly2[0].vertex, newpoly2[1].vertex)
                                newpoly2[len(newpoly2)-1].normal = get_normal(newpoly2[len(newpoly2)-2].vertex, newpoly2[len(newpoly2)-1].vertex, newpoly2[0].vertex)
                                newpoly1[0].normal = get_normal(newpoly1[len(newpoly1)-1].vertex, newpoly1[0].vertex, newpoly1[1].vertex)
                                newpoly1[len(newpoly1)-1].normal = get_normal(newpoly1[len(newpoly1)-2].vertex, newpoly1[len(newpoly1)-1].vertex, newpoly1[0].vertex)
                                newloopvertices.append(newpoly2)
                                # renumber additions
                                newadditions = {}
                                for keya1, va1 in additions.items():
                                    a1 = fromlistkey(keya1)
                                    newa1 = [a1[0], a1[1]]
                                    newva1 = [va1[0], va1[1], va1[2]]
                                    if a1[0] == a[0] and a1[1] == a[1]:
                                        pass
                                    elif a1[0] == a[0]:
                                        if newva1[1] in newpoly2_indexes:
                                            newva1[1] = newpoly2_indexes[newva1[1]]
                                            newva1[0] = len(newloopvertices)-1
                                        if newa1[1] in newpoly2_indexes:
                                            newa1[1] = newpoly2_indexes[newa1[1]]
                                            newa1[0] = len(newloopvertices)-1
                                        if newva1[1] in newpoly1_indexes:
                                            newva1[1] = newpoly1_indexes[newva1[1]]
                                        if newa1[1] in newpoly1_indexes:
                                            newa1[1] = newpoly1_indexes[newa1[1]]
                                        newadditions[listkey(newa1)] = newva1
                                    else:
                                        newadditions[listkey(newa1)] = newva1
                                print("old additions: " + str(additions))
                                additions = newadditions
                                newloopvertices[a[0]] = newpoly1
                                print("Broke loop into two: " + str(newpoly1) + " and " + str(newpoly2))
                                print("new additions: " + str(newadditions))
                                break
                            else:
                                # two different loops intersect
                                # create a single loop preserving direction:  (va[0], va[1]...va[1])...(a[0]...next(a[1])...a[1])
                                # this case should never occur
                                print("WARNING: joining loops should never occur!\n")
                                polya_indexes = {}
                                polyva_indexes = {}
                                newpoly = []
                                recalc = []
                                i = va[1]
                                polyva_indexes[va[1]] = len(newpoly)
                                recalc.append(len(newpoly))
                                newpoly.append(newloopvertices[va[0]][va[1]])
                                i = i + 1
                                if i == len(newloopvertices[va[0]]):
                                    i = 0
                                while i != va[1]:
                                    polyva_indexes[i] = len(newpoly)
                                    newpoly.append(newloopvertices[va[0]][i])
                                    i = i + 1
                                    if i == len(newloopvertices[va[0]]):
                                        i = 0
                                newvert = bm.verts.new(newloopvertices[va[0]][va[1]].vertex);
                                newVertexNormal = VertexNormal(newvert.co, (0.0, 0.0, 0.0), newvert)
                                recalc.append(len(newpoly))
                                newpoly.append(newVertexNormal)
                                i = a[1] + 1
                                if i == len(newloopvertices[a[0]]):
                                    i = 0
                                polya_indexes[i] = len(newpoly)
                                recalc.append(len(newpoly))
                                newpoly.append(newloopvertices[a[0]][i])
                                si = i
                                i = i + 1
                                if i == len(newloopvertices[a[0]]):
                                    i = 0
                                while i != si:
                                    polyva_indexes[i] = len(newpoly)
                                    newpoly.append(newloopvertices[a[0]][i])
                                    i = i + 1
                                    if i == len(newloopvertices[a[0]]):
                                        i = 0
                                recalc.append(len(newpoly)-1)
                                for r in recalc:
                                    pr = r - 1
                                    if pr < 0:
                                        pr = len(newpoly)-1
                                    nr = r + 1
                                    if nr == len(newpoly):
                                        nr = 0
                                    newpoly[r].normal = get_normal(newpoly[pr].vertex, newpoly[r].vertex, newpoly[nr].vertex)
                                # renumber additions
                                newadditions = {}
                                newloopindex = a[0]
                                removedloopindex = va[0]
                                if va[0] < newloopindex:
                                    newloopindex = va[0]
                                    removedloopindex = a[0]
                                for keya1, va1 in additions.items():
                                    a1 = fromlistkey(keya1)
                                    newa1 = [a1[0], a1[1]]
                                    newva1 = [va1[0], va1[1], va1[2]]
                                    if a1[0] == a[0] and a1[1] == a[1]:
                                        pass
                                    else:
                                        if newva1[0] == va[0]:
                                            newva1[1] = polyva_indexes[newva1[1]]
                                            newva1[0] = newloopindex
                                        if a1[0] == a[0]:
                                            newa1[1] = polya_indexes[newa1[1]]
                                            newa1[0] = newloopindex
                                        if newa1[0] > removedloopindex:
                                            newa1[0] = newa1[0] - 1
                                        if newva1[0] > removedloopindex:
                                            newva1[0] = newva1[0] - 1
                                        newadditions[listkey([newa1])] = newva1
                                print("old additions: " + str(additions))
                                additions = newadditions
                                newloopvertices[newloopindex] = newpoly
                                del newloopvertices[removedloopindex]
                                print("Joined loop into one: " + str(newpoly))
                                print("new additions: " + str(newadditions))
                                break
                else:    
                    # also remove overlapping points if their distance isn't increasing
                    for polyindex in range(0, len(newloopvertices)):
                        poly = newloopvertices[polyindex]
                        for i in range(0, len(poly)):
                            if not(listkey([polyindex,i]) in removals):
                                ni = i + 1
                                if ni >= len(poly):
                                    ni = 0
                                if i != ni and not(listkey([polyindex,ni]) in removals):
                                    d = vmin(poly[i].vertexObject.co, poly[ni].vertexObject.co)
                                    if vnorm2(d) < 1e-10:
                                        removals[listkey([polyindex,i])] = [6,0]
                    # process removals
                    orig_removals = {}
                    for k,v in removals.items():
                        orig_removals[k] = v
                    for polyindex in range(0, len(newloopvertices)):
                        newpoly = newloopvertices[polyindex]
                        faces = loopfaces[polyindex]
                        i = 0
                        while i < len(newpoly):
                            if listkey([polyindex, i]) in removals:
                                ni = i + 1
                                if ni == len(newpoly):
                                    ni = 0
                                # recalculate the normal of the neighbors (previous and next)
                                pi = i - 1
                                if pi < 0:
                                    pi = len(newpoly)-1
                                nni = ni+1
                                if nni == len(newpoly):
                                    nni = 0
                                ppi = pi - 1
                                if ppi < 0:
                                    ppi = len(newpoly)-1
                                pn = get_normal(newpoly[ppi].vertexObject.co, newpoly[pi].vertexObject.co, newpoly[ni].vertexObject.co)
                                nn = get_normal(newpoly[pi].vertexObject.co, newpoly[ni].vertexObject.co, newpoly[nni].vertexObject.co)
                                # do not move away fixed points
                                if vnorm2(newpoly[pi].normal) != 0.0:
                                    newpoly[pi].normal = pn
                                if vnorm2(newpoly[ni].normal) != 0.0:
                                    newpoly[ni].normal = nn
                                bm.verts.remove(newpoly[i].vertexObject)
                                del newpoly[i]
                                newni = i
                                if newni == len(newpoly):
                                    # when i was the last
                                    newni = 0
                                print("Processing removal: " + str([polyindex, i]) + ", ni = " + str(ni) + ", pi = " + str(pi) + ", ppi = " + str(ppi) + ", nni = " + str(nni) + ", newni = " + str(newni))
                                newfaces = []
                                for fi in range(0, len(faces)):
                                    if (faces[fi][2] == 1 and faces[fi][3] == ni and faces[fi][4] == 1 and faces[fi][5] == i) or (pi == i) or (ni == i):
                                        print("Loop " + str(polyindex) + ", Face removal: " + str(faces[fi]))           
                                        pass
                                    else:
                                        print("Loop " + str(polyindex) + ", old face: " + str(faces[fi]))
                                        newface = [faces[fi][0],faces[fi][1], faces[fi][2], faces[fi][3], faces[fi][4], faces[fi][5]]
                                        if faces[fi][2] == 1 and faces[fi][3] > i:
                                            newface[3] = newface[3] -1
                                        if faces[fi][4] == 1 and faces[fi][5] > i:
                                            newface[5] = newface[5] -1
                                        if faces[fi][2] == 1 and faces[fi][3] == i:
                                            newface[3] = newni
                                        if faces[fi][4] == 1 and faces[fi][5] == i:
                                            newface[5] = newni
                                        print("Loop " + str(polyindex) + ", new face: " + str(newface))
                                        if faces[fi][2] == 1 and newface[3] == newface[5]:
                                            print("So many removals that face became degenerate: " + str(faces[fi]))
                                            pass
                                        else:
                                            newfaces.append(newface)
                                faces = newfaces
                                loopfaces[polyindex] = faces
                                # update removals
                                newremovals = {}
                                for rl,rv in removals.items():
                                    r = fromlistkey(rl)
                                    newr = [r[0], r[1]]
                                    if r[0] == polyindex:
                                        if r[1] == i:
                                            pass
                                        else:
                                            if r[1] > i:
                                                newr[1] = newr[1] - 1
                                            newremovals[listkey(newr)] = rv
                                    else:
                                        newremovals[listkey(newr)] = rv
                                removals = newremovals
                                # update closing faces
                                # newclosingfaces = {}
                                # for fk, f in closingfaces.items():
                                #     newclosingface = []
                                #     for fv in f:
                                #         if fv[0] == polyindex and fv[1] == i:
                                #             newni = ni
                                #             if newni > i:
                                #                 newni = newni - 1
                                #             if not([fv[0], newni] in newclosingface):
                                #                 newclosingface.append([fv[0], newni])
                                #         elif fv[0] == polyindex and fv[1] > i:
                                #             if not([fv[0], fv[1]-1] in newclosingface):
                                #                 newclosingface.append([fv[0], fv[1]-1])
                                #         else:
                                #             newclosingface.append([fv[0], fv[1]])
                                #     if len(newclosingface) > 2:
                                #         newclosingfaces[str(newclosingface)] = newclosingface
                                # closingfaces = newclosingfaces                            
                            else:
                                i = i+1
                    # detect degenerate cases:
                    # dotprod(pv' - v' , nv' - v') > 0, crossprod(pv' - v', nv' - v') = 0   -->   remove v'
                    #
                    removals = {}
                    for newpolyindex in range(0, len(newloopvertices)): 
                        newpoly = newloopvertices[newpolyindex]
                        i = 0
                        while i < len(newpoly):
                            ni = i + 1
                            if ni == len(newpoly):
                                ni = 0
                                # recalculate the normal of the neighbors (previous and next)
                            pi = i - 1
                            if pi < 0:
                                pi = len(newpoly)-1
                            pvv = vmin(newpoly[pi].vertexObject.co, newpoly[i].vertexObject.co)
                            nvv = vmin(newpoly[ni].vertexObject.co, newpoly[i].vertexObject.co)
                            if vdot(pvv, nvv) > 0 and vnorm2(vcross(pvv, nvv)) < 1e-5 * vnorm2(pvv) * vnorm2(nvv) and not(listkey([newpolyindex,pi]) in removals) and not(listkey([newpolyindex,ni]) in removals) :
                                removals[listkey([newpolyindex, i])] = [5, 0]
                            i = i + 1
                    # process removals again
                    for polyindex in range(0, len(newloopvertices)):
                        newpoly = newloopvertices[polyindex]
                        faces = loopfaces[polyindex]
                        i = 0
                        while i < len(newpoly):
                            if listkey([polyindex, i]) in removals:
                                ni = i + 1
                                if ni == len(newpoly):
                                    ni = 0
                                # recalculate the normal of the neighbors (previous and next)
                                pi = i - 1
                                if pi < 0:
                                    pi = len(newpoly)-1
                                nni = ni+1
                                if nni == len(newpoly):
                                    nni = 0
                                ppi = pi - 1
                                if ppi < 0:
                                    ppi = len(newpoly)-1
                                pn = get_normal(newpoly[ppi].vertexObject.co, newpoly[pi].vertexObject.co, newpoly[ni].vertexObject.co)
                                nn = get_normal(newpoly[pi].vertexObject.co, newpoly[ni].vertexObject.co, newpoly[nni].vertexObject.co)
                                # do not move away fixed points
                                if vnorm2(newpoly[pi].normal) != 0.0:
                                    newpoly[pi].normal = pn
                                if vnorm2(newpoly[ni].normal) != 0.0:
                                    newpoly[ni].normal = nn
                                bm.verts.remove(newpoly[i].vertexObject)
                                del newpoly[i]
                                newni = i
                                if newni == len(newpoly):
                                    # when i was the last
                                    newni = 0
                                print("Processing removal: " + str([polyindex, i]) + ", ni = " + str(ni) + ", pi = " + str(pi) + ", ppi = " + str(ppi) + ", nni = " + str(nni) + ", newni = " + str(newni))
                                newfaces = []
                                for fi in range(0, len(faces)):
                                    if (faces[fi][2] == 1 and faces[fi][3] == ni and faces[fi][4] == 1 and faces[fi][5] == i)  or (pi == i) or (ni == i):
                                        print("Loop " + str(polyindex) + ", Face removal: " + str(faces[fi]))           
                                        pass
                                    else:
                                        print("Loop " + str(polyindex) + ", old face: " + str(faces[fi]))
                                        newface = [faces[fi][0],faces[fi][1], faces[fi][2], faces[fi][3], faces[fi][4], faces[fi][5]]
                                        if faces[fi][2] == 1 and faces[fi][3] > i:
                                            newface[3] = newface[3] -1
                                        if faces[fi][4] == 1 and faces[fi][5] > i:
                                            newface[5] = newface[5] -1
                                        if faces[fi][2] == 1 and faces[fi][3] == i:
                                            newface[3] = newni
                                        if faces[fi][4] == 1 and faces[fi][5] == i:
                                            newface[5] = newni
                                        print("Loop " + str(polyindex) + ", new face: " + str(newface))
                                        if faces[fi][2] == 1 and newface[3] == newface[5]:
                                            print("So many removals that face became degenerate: " + str(faces[fi]))
                                            pass
                                        else:
                                            newfaces.append(newface)
                                faces = newfaces
                                loopfaces[polyindex] = faces
                                # update removals
                                newremovals = {}
                                for rl,rv in removals.items():
                                    r = fromlistkey(rl)
                                    newr = [r[0], r[1]]
                                    if r[0] == polyindex:
                                        if r[1] == i:
                                            pass
                                        else:
                                            if r[1] > i:
                                                newr[1] = newr[1] - 1
                                            newremovals[listkey(newr)] = rv
                                    else:
                                        newremovals[listkey(newr)] = rv
                                removals = newremovals
                                # update closing faces
                                # newclosingfaces = {}
                                # for fk, f in closingfaces.items():
                                #     newclosingface = []
                                #     for fv in f:
                                #         if fv[0] == polyindex and fv[1] == i:
                                #             newni = ni
                                #             if newni > i:
                                #                 newni = newni - 1
                                #             if not([fv[0], newni] in newclosingface):
                                #                 newclosingface.append([fv[0], newni])
                                #         elif fv[0] == polyindex and fv[1] > i:
                                #             if not([fv[0], fv[1]-1] in newclosingface):
                                #                 newclosingface.append([fv[0], fv[1]-1])
                                #         else:
                                #             newclosingface.append([fv[0], fv[1]])
                                #     if len(newclosingface) > 2:
                                #         newclosingfaces[str(newclosingface)] = newclosingface
                                # closingfaces = newclosingfaces                            
                            else:
                                i = i+1                    
                    # add faces
                    for polyindex in range(0, len(loopvertices)): 
                        faces = loopfaces[polyindex]
                        poly = loopvertices[polyindex]
                        newpoly = newloopvertices[polyindex]
                        for face in faces:
                            dprint("Loop " + str(polyindex) + ", Face creation: " + str(face))                            
                            if face[2] == 0:
                                if face[1] >= len(poly) or face[3] >= len(poly) or face[5] >= len(newpoly):
                                    print("Incorrect face index: len(poly) = " + str(len(poly)) + ", len(newpoly) = " + str(len(newpoly)) + ", face = " + str(face))
                                else:
                                    addface( bm, (poly[face[1]].vertexObject, poly[face[3]].vertexObject, newpoly[face[5]].vertexObject), sampleFace)
                            else:
                                if face[3] != face[5]:
                                    addface( bm, (poly[face[1]].vertexObject, newpoly[face[3]].vertexObject, newpoly[face[5]].vertexObject), sampleFace)
                # for all newly added vertices check the vertex connected at the below level in the zdepth+normal direction
                # for each check the triangles it is connected to
                #    if moving that vertex to the location of the newly added vertex wouldn't change the normal of any face, then:
                #    delete all faces which are connected both to the newly added vertex and the one below
                #    recreate all remaining faces of the below vertex with the new vertex (delete the old ones)
                #    delete the below vertex
                for polyindex in range(0, len(newloopvertices)): 
                    newpoly = newloopvertices[polyindex]
                    for i in range(0, len(newpoly)):
                        dprint("Checking loop " + str(polyindex) + ", vertex " + str(i) + " if faces can be merged: " + str(newpoly[i].vertexObject.co))

                        zdepth_plus_normal = vadd(zdepth, newpoly[i].normal)
                        dprint("zdepth_plus_normal = " + str(zdepth_plus_normal) + ", normalized: " + str(vscale(zdepth_plus_normal, 1.0/math.sqrt(vnorm2(zdepth_plus_normal)))))
                        below = []
                        zbelow = []
                        for fk, f in enumerate(newpoly[i].vertexObject.link_faces):
                            for vk, v in enumerate(f.verts):
                                diff = vmin(v.co,newpoly[i].vertexObject.co)
                                if vnorm2(diff) > 0:
                                    dprint("diff = " + str(diff) + ", normalized = " + str(vscale(diff, 1.0/math.sqrt(vnorm2(diff)))))
                                else:
                                    dprint("diff = " + str(diff))
                                if vnorm2(diff) > 0.0:
                                    if vnorm2(newpoly[i].normal) > 0.0 and 1. - vdot(diff, zdepth_plus_normal)*vdot(diff, zdepth_plus_normal) / (vdot(diff, diff)*vdot(zdepth_plus_normal, zdepth_plus_normal)) < 1e-5:
                                        if not v in below:
                                            below.append(v)
                                        dprint("below: " + str(v.co))
                                    # look for below vectors in the zdepth direction too, after the vertex shift stopped
                                    elif 1. - vdot(diff, zdepth)*vdot(diff, zdepth) / (vdot(diff, diff)*vdot(zdepth, zdepth)) < 1e-5:
                                        if not v in zbelow:
                                            zbelow.append(v)
                                        dprint("zbelow: " + str(v.co))
                                    else:
                                        dprint("not below: " + str(v.co) + ", by: " + str(1. - vdot(diff, zdepth_plus_normal)*vdot(diff, zdepth_plus_normal) / (vdot(diff, diff)*vdot(zdepth_plus_normal, zdepth_plus_normal))))
                        if len(below) + len(zbelow) > 0:
                            dprint("Found the vector(s) below: " + str(below) + ", zbelow: " + str(zbelow))
                            ok = True
                            for bv in below:
                                for fk, f in enumerate(bv.link_faces):
                                    # if the normal of the face is perpendicular to zdepth_plus_normal then we are fine to move 'below' to 'top'
                                    # normal_update is called on adding each new face in addface, so they have a normal.
                                    if vnorm2(f.normal) > 0 and vdot(f.normal, zdepth_plus_normal)*vdot(f.normal, zdepth_plus_normal) / vnorm2(f.normal)/vnorm2(zdepth_plus_normal) > 1e-3:
                                        ok = False
                                        dprint("not perpendicular by: " + str(vdot(f.normal, zdepth_plus_normal)*vdot(f.normal, zdepth_plus_normal) / vnorm2(f.normal)/vnorm2(zdepth_plus_normal)))
                            for bv in zbelow:
                                for fk, f in enumerate(bv.link_faces):
                                    # if the normal of the face is perpendicular to zdepth then we are fine to move 'below' to 'top'
                                    # normal_update is called on adding each new face in addface, so they have a normal.
                                    if vnorm2(f.normal) > 0 and vdot(f.normal, zdepth)*vdot(f.normal, zdepth) / vnorm2(f.normal)/vnorm2(zdepth) > 1e-3:
                                        ok = False
                                        dprint("not perpendicular by: " + str(vdot(f.normal, zdepth)*vdot(f.normal, zdepth) / vnorm2(f.normal)/vnorm2(zdepth)))
                            dprint("If faces are perpendicular: " + str(ok))
                            if ok:
                                toRemove = []
                                toCreate = []
                                for bv in below + zbelow:
                                    for fk, f in enumerate(bv.link_faces):
                                        remove = False
                                        for vk, v in enumerate(f.verts):
                                            d = vmin(v.co, newpoly[i].vertexObject.co)
                                            if vnorm2(d) < 1e-15:
                                                dprint("Found face with both below and new: " + str(f.verts))
                                                toRemove.append(f)
                                                remove = True
                                                break
                                        if not remove:
                                            toRemove.append(f)
                                            fverts = []
                                            for vk, v in enumerate(f.verts):
                                                d = vmin(v.co, bv.co)
                                                if vnorm2(d) < 1e-15:
                                                    vertToAdd = newpoly[i].vertexObject
                                                else:
                                                    vertToAdd = v
                                                if not (vertToAdd in fverts):
                                                    fverts.append(vertToAdd)
                                            if len(fverts) > 2:
                                                toCreate.append(fverts)
                                            dprint("Face recreated, was: " + str(f.verts) + ", will be: " + str(fverts))
                                for f in toRemove:
                                    if f.is_valid:
                                        bm.faces.remove(f)
                                for verts in toCreate:
                                    addface(bm, verts, sampleFace)
                                if bv.is_valid:
                                    bm.verts.remove(bv)
                loopvertices = []
                for polyindex in range(0, len(newloopvertices)):
                    poly = newloopvertices[polyindex]
                    if len(poly) > 2:
                        loopvertices.append(poly)
                if len(loopvertices) == 0:
                    break
                # do limited loops only
                print("Shift loop completed: " + str(iter))
                with open('/tmp/loop-' + str(iter) + '.txt', 'w') as f:
                    f.write("# original additions: " + str(orig_additions) + "\n")
                    f.write("# original removals: " + str(orig_removals) + "\n")
                    f.write("# final additions: " + str(additions) + "\n")
                    f.write("# final removals: " + str(removals) + "\n")
                    for polyindex in range(0, len(loopvertices)): 
                        poly = loopvertices[polyindex]
                        for i in range(0, len(poly)):
                            ni = i + 1
                            if ni == len(poly):
                                ni = 0
                            length = math.sqrt(vnorm2(vmin(poly[ni].vertex, poly[i].vertex)))
                            f.write('{0:.8g} {1:.8g} {2:.8g} {3:.8g} {4:.8g} {5:.8g} {6:.8g} {7:.8g} {8:.8g} {9} \n'.format(poly[i].vertex[0], poly[i].vertex[1], poly[i].vertex[2], poly[ni].vertex[0] - poly[i].vertex[0], poly[ni].vertex[1] - poly[i].vertex[1], poly[ni].vertex[2] - poly[i].vertex[2], length*poly[i].normal[0],  length*poly[i].normal[1], length*poly[i].normal[2], str(polyindex) + '-' + str(i)))
                iter = iter + 1
                if iter > 2000:
                    break
            #            ret = bmesh.ops.dissolve_limit(bm, angle_limit=1e-5, use_dissolve_boundaries=False, verts=bm.verts, edges=bm.edges, delimit=1)
            print("Shift loops completed.")
            # add closing faces
            # for fk, f in closingfaces.items():
            #     vertices = []
            #     for fv in f:
            #         vertices.append(loopvertices[fv[0]][fv[1]].vertexObject)
            #     newface = bm.faces.new(vertices, sampleFace)
            #     newface.normal_update()
            #     if vdot(newface.normal, zdepth) > 0:
            #         newface.normal_flip()
            bmesh.update_edit_mesh(img)
