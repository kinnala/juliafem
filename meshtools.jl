#
# Some mesh tools for finite elements.
#
# Author: Tom Gustafsson, 11.2.2014
# Licensed under GPLv3
#

# Allows to use the plotting commands
# of matplotlib
using PyPlot

type Mesh
    # Mesh type for storing
    # triangular two-dimesional
    # finite element meshes

    p::Array{Float64,2}
    t::Array{Int64,2}
    edges::Array{Int64,2}
    t2e::Array{Int64,2}

    function Mesh(p,t)
        # Constructor to build the edge
        # numbering etc
        # TODO Add the edge2triangle mapping wizardry

        t = sort(t,1)
        edges = Array(Int64,2,3*size(t,2))
        t2e = Array(Int64,3,size(t,2))
        check = Dict{(Int64,Int64),Int64}()
        N = 1;
        for k=1:size(t,2)
            n1=t[1,k]
            n2=t[2,k]
            n3=t[3,k]
            t2e1 = 0
            t2e2 = 0
            t2e3 = 0
            
            if !haskey(check,(n1,n2))
                edges[1,N]=n1;
                edges[2,N]=n2;
                t2e1 = N;
                check[(n1,n2)]=N
                N += 1;
            else
                t2e1 = check[(n1,n2)]
            end
            if !haskey(check,(n2,n3))
                edges[1,N]=n2;
                edges[2,N]=n3;
                t2e2 = N;
                check[(n2,n3)]=N
                N += 1;
            else
                t2e2 = check[(n2,n3)]
            end
            if !haskey(check,(n1,n3))
                edges[1,N]=n1;
                edges[2,N]=n3;
                t2e3 = N;
                check[(n1,n3)]=N
                N += 1;
            else
                t2e3 = check[(n1,n3)]
            end
            t2e[1,k] = t2e1
            t2e[2,k] = t2e2
            t2e[3,k] = t2e3
        end
        new(p,t,edges[:,1:(N-1)],t2e)
    end
end

# Default constructor
Mesh() = Mesh([0 0; 1 0; 0 1; 1 1]',[1 2 3; 2 3 4]')

function trimesh(mesh)
    # Plot a finite element mesh
    figure()
    tx = Array(Float64,size(mesh.t,1)+1,size(mesh.t,2))
    ty = Array(Float64,size(mesh.t,1)+1,size(mesh.t,2))
    for k=1:size(mesh.t,2)
        tx[:,k] = mesh.p'[mesh.t[[1,2,3,1],k],1]
        ty[:,k] = mesh.p'[mesh.t[[1,2,3,1],k],2]
    end
    plot(tx,ty,"k-")
end

function triplot(mesh,u,N=10)
    # Draw a contour plot of "u" on "mesh"
    figure()
    tricontour(mesh.p'[:,1],mesh.p'[:,2],mesh.t'-1,u,N)
end

function trisurf(mesh,u)
    figure()
    # TODO figure a way to add the mesh and some
    # colormap here.
    plot_trisurf(mesh.p'[:,1],mesh.p'[:,2],u,cmap=ColorMap("jet"))
end

function refinetri(mesh::Mesh)
    Np = size(mesh.p,2)
    Ne = size(mesh.edges,2)
    Nt = size(mesh.t,2)
    # New points
    rp = Array(Float64,2,Np+Ne)
    for i=1:Np
        rp[1,i] = mesh.p[1,i]
        rp[2,i] = mesh.p[2,i]
    end
    for j=1:Ne
        rp[1,Np+j] = 0.5*(mesh.p[1,mesh.edges[1,j]]+mesh.p[1,mesh.edges[2,j]])
        rp[2,Np+j] = 0.5*(mesh.p[2,mesh.edges[1,j]]+mesh.p[2,mesh.edges[2,j]])
    end
    # New triangles
    rt = Array(Int64,3,4*Nt)
    for k=1:Nt
        rt[1,k] = mesh.t[1,k]
        rt[1,k+Nt] = mesh.t[2,k]
        rt[1,k+2*Nt] = mesh.t[3,k]
        rt[1,k+3*Nt] = mesh.t2e[1,k]+Np
        rt[2,k] = mesh.t2e[1,k]+Np;
        rt[2,k+Nt] = mesh.t2e[1,k]+Np;
        rt[2,k+2*Nt] = mesh.t2e[3,k]+Np;
        rt[2,k+3*Nt] = mesh.t2e[2,k]+Np;
        rt[3,k] = mesh.t2e[3,k]+Np;
        rt[3,k+Nt] = mesh.t2e[2,k]+Np;
        rt[3,k+2*Nt] = mesh.t2e[2,k]+Np;
        rt[3,k+3*Nt] = mesh.t2e[3,k]+Np;
    end

    return Mesh(rp,rt)
end

function meshtest()
    mesh = Mesh()
    for i=1:4
        mesh = refinetri(mesh)
    end
    f(x,y) = x.*y
    triplot(mesh,f(mesh.p'[:,1],mesh.p'[:,2]),20)
end
