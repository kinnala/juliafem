module Fem
#
# Simple finite element assembly
# using piecewise linear triangular
# elements.
#
# Author: Tom Gustafsson, 30.3.2014
# Licensed under GPLv3
#

include("meshtools.jl")

macro bilin_asm(bilin_form,mesh)
    quote
        ##
        ## Assembly of a bilinear form using P1 elements.
        ## Results in a sparse matrix K_ij.
        ##

        # Store triplets (i,j,K_ij)
        I = Array(Int64,0)
        J = Array(Int64,0)
        V = Array(Float64,0)

        # Local basis functions and their gradients
        phi_hat = [(x,y)->1-x-y, (x,y)->x, (x,y)->y]
        dphi_hat = [-1 -1; 1 0; 0 1]'

        # Dunavant quadrature rule. This is in fact
        # good for 2nd order polynomials on triangle.
        qp = [0.5 0; 0.5 0.5; 0 0.5]'
        qw = [1/6,1/6,1/6]
        Nqp = size(qp,2)

        # Array for storing the pre-evaluated basis functions
        # at the quadrature points
        phi_hat_qp = Array(Float64,3,Nqp)

        # Pre-evaluate basis functions at (local) QP's
        for q=1:Nqp
            for i=1:3
              phi_hat_qp[i,q] = phi_hat[i](qp[1,q],qp[2,q])
            end
        end

        # Allocate temporary storage for
        #   A         - affine mapping
        #   invAt     - inverse-transpose of affine mapping
        #   n1,n2,n3  - node locations
        #   gqp       - global quadrature points
        #   dphi      - global basis function values
        #   absdetA   - absolute value of determinant of affine map
        #   detA      - determinant of affine map
        A = Array(Float64,2,2)
        invAt = Array(Float64,2,2)
        n1 = Array(Float64,2)
        n2 = Array(Float64,2)
        n3 = Array(Float64,2)
        gqp = Array(Float64,2,Nqp)
        dphi = Array(Float64,size(dphi_hat,1),size(dphi_hat,2))
        absdetA::Float64 = -1.0
        detA::Float64 = -1.0

        for k=1:size($mesh.t,2)
            # Find global node locations
            n1[1] = $mesh.p[1,$mesh.t[1,k]]
            n1[2] = $mesh.p[2,$mesh.t[1,k]]
            n2[1] = $mesh.p[1,$mesh.t[2,k]]
            n2[2] = $mesh.p[2,$mesh.t[2,k]]
            n3[1] = $mesh.p[1,$mesh.t[3,k]]
            n3[2] = $mesh.p[2,$mesh.t[3,k]]
            # Build the affine mapping
            A[1,1] = n2[1]-n1[1]
            A[1,2] = n3[1]-n1[1]
            A[2,1] = n2[2]-n1[2]
            A[2,2] = n3[2]-n1[2]
            # Calculate global quadrature points
            for l=1:Nqp
                gqp[1,l] = A[1,1]*qp[1,l]+A[1,2]*qp[2,l]+n1[1]
                gqp[2,l] = A[2,1]*qp[1,l]+A[2,2]*qp[2,l]+n1[2]
            end
            # The determinant of A and its absolute value
            detA = A[1,1]*A[2,2]-A[1,2]*A[2,1]
            if detA<0
                absdetA = -detA
            else
                absdetA = detA
            end
            # Inverse of A
            invAt[1,1] = A[2,2]/detA
            invAt[2,1] = -A[1,2]/detA
            invAt[1,2] = -A[2,1]/detA
            invAt[2,2] = A[1,1]/detA
            # Calculate global gradients using invAt
            for l=1:size(dphi_hat,2)
                dphi[1,l] = invAt[1,1]*dphi_hat[1,l]+invAt[1,2]*dphi_hat[2,l]
                dphi[2,l] = invAt[2,1]*dphi_hat[1,l]+invAt[2,2]*dphi_hat[2,l]
            end

            # Loop over the local stiffness matrix
            for j=1:3
                ux::Float64 = dphi[1,j]
                uy::Float64 = dphi[2,j]
                jx::Int64 = $mesh.t[j,k]
                for i=1:3
                    vx::Float64 = dphi[1,i]
                    vy::Float64 = dphi[2,i]
                    ix::Int64 = $mesh.t[i,k]
                    # Loop over quadrature points
                    for q=1:Nqp
                        x = gqp[1,q]
                        y = gqp[2,q]
                        u = phi_hat_qp[j,q]
                        v = phi_hat_qp[i,q]
                        # Local stiffness matrix values are
                        # stored in V and the corresponding
                        # global indices in I and J.
                        # Such lists can be directly feeded
                        # to the sparse matrix constructor.
                        push!(I,ix)
                        push!(J,jx)
                        push!(V,($bilin_form)*qw[q]*absdetA)
                    end
                end
            end
        end
        # Finally build the stiffness matrix
        sparse(I,J,V)
    end
end

macro bilin_asm_mat(bilin_form,mesh,mat)
    quote
        ##
        ## Assembly of a (non-linear) bilinear form using P1 elements.
        ## Results in a sparse matrix K_ij.
        ## Additional parameter "nl" is a vector of nodal values
        ## of some parameter, which can be included in the bilinear
        ## form by using the variable "a".
        ##

        # Store triplets (i,j,K_ij)
        I = Array(Int64,0)
        J = Array(Int64,0)
        V = Array(Float64,0)

        # Local basis functions and their gradients
        phi_hat = [(x,y)->1-x-y, (x,y)->x, (x,y)->y]
        dphi_hat = [-1.0 -1.0; 1.0 0.0; 0.0 1.0]'

        # Dunavant quadrature rule. This is in fact
        # good for 2nd order polynomials on triangle.
        qp = [0.5 0; 0.5 0.5; 0 0.5]'
        qw = [1/6,1/6,1/6]
        Nqp = size(qp,2)

        # Array for storing the pre-evaluated basis functions
        # at the quadrature points
        phi_hat_qp = Array(Float64,3,Nqp)

        # Pre-evaluate basis functions at (local) QP's
        for q=1:Nqp
            for i=1:3
              phi_hat_qp[i,q] = phi_hat[i](qp[1,q],qp[2,q])
            end
        end

        # Allocate temporary storage for
        #   A           - affine mapping
        #   invAt       - inverse-transpose of affine mapping
        #   n1,n2,n3    - node locations
        #   n1i,n2i,n3i - node indices
        #   gqp         - global quadrature points
        #   dphi        - global basis function values
        #   absdetA     - absolute value of determinant of affine map
        #   detA        - determinant of affine map
        A = Array(Float64,2,2)
        invAt = Array(Float64,2,2)
        n1 = Array(Float64,2)
        n2 = Array(Float64,2)
        n3 = Array(Float64,2)
        n1i::Int64 = -1
        n2i::Int64 = -1
        n3i::Int64 = -1
        gqp = Array(Float64,2,Nqp)
        dphi = Array(Float64,size(dphi_hat,1),size(dphi_hat,2))
        absdetA::Float64 = -1.0
        detA::Float64 = -1.0

        nlcalc = $mat

        for k=1:size($mesh.t,2)
            # Find global node indices
            n1i = $mesh.t[1,k]
            n2i = $mesh.t[2,k]
            n3i = $mesh.t[3,k]
            # Find the parameter values
            a1::Float64 = nlcalc[n1i,1]
            a2::Float64 = nlcalc[n2i,1]
            a3::Float64 = nlcalc[n3i,1]
            # Find global node locations
            n1[1] = $mesh.p[1,n1i]
            n1[2] = $mesh.p[2,n1i]
            n2[1] = $mesh.p[1,n2i]
            n2[2] = $mesh.p[2,n2i]
            n3[1] = $mesh.p[1,n3i]
            n3[2] = $mesh.p[2,n3i]
            # Build the affine mapping
            A[1,1] = n2[1]-n1[1]
            A[1,2] = n3[1]-n1[1]
            A[2,1] = n2[2]-n1[2]
            A[2,2] = n3[2]-n1[2]
            # Calculate global quadrature points
            for l=1:Nqp
                gqp[1,l] = A[1,1]*qp[1,l]+A[1,2]*qp[2,l]+n1[1]
                gqp[2,l] = A[2,1]*qp[1,l]+A[2,2]*qp[2,l]+n1[2]
            end
            # The determinant of A and its absolute value
            detA = A[1,1]*A[2,2]-A[1,2]*A[2,1]
            if detA<0
                absdetA = -detA
            else
                absdetA = detA
            end
            # Inverse of A
            invAt[1,1] = A[2,2]/detA
            invAt[2,1] = -A[1,2]/detA
            invAt[1,2] = -A[2,1]/detA
            invAt[2,2] = A[1,1]/detA
            # Calculate global gradients using invAt
            for l=1:size(dphi_hat,2)
                dphi[1,l] = invAt[1,1]*dphi_hat[1,l]+invAt[1,2]*dphi_hat[2,l]
                dphi[2,l] = invAt[2,1]*dphi_hat[1,l]+invAt[2,2]*dphi_hat[2,l]
            end

            # Loop over the local stiffness matrix
            for j=1:3
                ux::Float64 = dphi[1,j]
                uy::Float64 = dphi[2,j]
                jx::Int64 = $mesh.t[j,k]
                for i=1:3
                    vx::Float64 = dphi[1,i]
                    vy::Float64 = dphi[2,i]
                    ix::Int64 = $mesh.t[i,k]
                    # Loop over quadrature points
                    for q=1:Nqp
                        x = gqp[1,q]
                        y = gqp[2,q]
                        u = phi_hat_qp[j,q]
                        v = phi_hat_qp[i,q]
                        a = a1+(a2-a1)*qp[1,q]+(a3-a1)*qp[2,q]
                        # Local stiffness matrix values are
                        # stored in V and the corresponding
                        # global indices in I and J.
                        # Such lists can be directly feeded
                        # to the sparse matrix constructor.
                        push!(I,ix)
                        push!(J,jx)
                        push!(V,($bilin_form)*qw[q]*absdetA)
                    end
                end
            end
        end
        # Finally build the stiffness matrix
        sparse(I,J,V)
    end
end

macro lin_asm(lin_form,mesh)
    quote
        ##
        ## Assembly of a linear form using P1 elements.
        ## Results in a vector f_i.
        ##
        I = Array(Int64,0)
        V = Array(Float64,0)

        # Local basis functions and their gradients
        phi_hat = [(x,y)->1-x-y, (x,y)->x, (x,y)->y]
        dphi_hat = [-1 -1; 1 0; 0 1]'

        # Dunavant quadrature rule. This is in fact
        # good for 2nd order polynomials on triangle.
        qp = [0.5 0; 0.5 0.5; 0 0.5]'
        qw = [1/6,1/6,1/6]
        Nqp = size(qp,2)

        # Array for storing the pre-evaluated basis functions
        # at the quadrature points
        phi_hat_qp = Array(Float64,3,Nqp)

        # Pre-evaluate basis functions at (local) QP's because
        # function evaluation is costly
        for q=1:Nqp
            for i=1:3
              phi_hat_qp[i,q] = phi_hat[i](qp[1,q],qp[2,q])
            end
        end

        # Allocate temporary storage for
        #   A         - affine mapping
        #   invAt     - inverse-transpose of affine mapping
        #   n1,n2,n3  - node locations
        #   gqp       - global quadrature points
        #   dphi      - global basis function values
        #   absdetA   - absolute value of determinant of affine map
        #   detA      - determinant of affine map
        A = Array(Float64,2,2)
        invAt = Array(Float64,2,2)
        n1 = Array(Float64,2)
        n2 = Array(Float64,2)
        n3 = Array(Float64,2)
        gqp = Array(Float64,2,Nqp)
        dphi = Array(Float64,size(dphi_hat,1),size(dphi_hat,2))
        absdetA::Float64 = -1.0
        detA::Float64 = -1.0

        for k=1:size($mesh.t,2)
            # Find nodes
            n1[1] = $mesh.p[1,$mesh.t[1,k]]
            n1[2] = $mesh.p[2,$mesh.t[1,k]]
            n2[1] = $mesh.p[1,$mesh.t[2,k]]
            n2[2] = $mesh.p[2,$mesh.t[2,k]]
            n3[1] = $mesh.p[1,$mesh.t[3,k]]
            n3[2] = $mesh.p[2,$mesh.t[3,k]]
            # Form the affine mapping using nodes
            A[1,1] = n2[1]-n1[1]
            A[1,2] = n3[1]-n1[1]
            A[2,1] = n2[2]-n1[2]
            A[2,2] = n3[2]-n1[2]
            # Calculate global quadrature points
            for l=1:Nqp
                gqp[1,l] = A[1,1]*qp[1,l]+A[1,2]*qp[2,l]+n1[1]
                gqp[2,l] = A[2,1]*qp[1,l]+A[2,2]*qp[2,l]+n1[2]
            end
            # The determinant of A and its absolute value
            detA = A[1,1]*A[2,2]-A[1,2]*A[2,1]
            if detA<0
                absdetA = -detA
            else
                absdetA = detA
            end
            # Inverse of A for transforming the gradients
            invAt[1,1] = A[2,2]/detA
            invAt[2,1] = -A[1,2]/detA
            invAt[1,2] = -A[2,1]/detA
            invAt[2,2] = A[1,1]/detA
            # Calculate global gradients
            for l=1:size(dphi_hat,2)
                dphi[1,l] = invAt[1,1]*dphi_hat[1,l]+invAt[1,2]*dphi_hat[2,l]
                dphi[2,l] = invAt[2,1]*dphi_hat[1,l]+invAt[2,2]*dphi_hat[2,l]
            end

            # Loop over the local load vector
            for i=1:3
                vx::Float64 = dphi[1,i]
                vy::Float64 = dphi[2,i]
                ix::Int64 = $mesh.t[i,k]
                # Loop over quadrature points
                for q=1:Nqp
                    x = gqp[1,q]
                    y = gqp[2,q]
                    v = phi_hat_qp[i,q]
                    # Local load vector values are
                    # stored in V and the corresponding
                    # global indices in I.
                    # Such lists can be directly feeded
                    # to the sparse matrix constructor.
                    push!(I,ix)
                    push!(V,($lin_form)*qw[q]*absdetA)
                end
            end
        end
        # Finally build the load vector
        sparsevec(I,V)
    end
end

function test_poisson()
    ##
    ## Test case solves the Poisson equation
    ##   u_xx + u_yy = 1,
    ## with the boundary condition
    ##   u = 0,
    ## on unit square.
    ##
    mesh = Mesh()
    for i=1:4
        mesh = refinetri(mesh)
    end
    S = @bilin_asm(ux*vx+uy*vy,mesh)
    f = @lin_asm(v,mesh)
    # Find Dirichlet node indices
    D = Array(Int64,0)
    for i=1:size(mesh.p,2)
        if mesh.p[1,i]==0 || mesh.p[1,i]==1 || mesh.p[2,i]==0 || mesh.p[2,i]==1
            push!(D,i)
        end
    end
    I = setdiff(1:size(mesh.p,2),D)
    # Solve the system
    u = Array(Float64,size(mesh.p,2))
    u = zeros(size(mesh.p,2))
    u[I]=cholfact(S[I,I])\f[I,1]

    # Check that the maximum value is correct
    exactmaxu = 0.07344576676891949
    eps = 0.00001
    if maximum(u)>exactmaxu+eps || maximum(u)<exactmaxu-eps
        return false
    end
    true
end

function runtests()
    print("test_poisson...")
    print(test_poisson())
    print("\n")
end

function play()
    mesh = Mesh()
    for i=1:4
        mesh = refinetri(mesh)
    end

    # Find Dirichlet node indices
    D = Array(Int64,0)
    for i=1:size(mesh.p,2)
        if mesh.p[1,i]==0 || mesh.p[1,i]==1 || mesh.p[2,i]==0 || mesh.p[2,i]==1
            push!(D,i)
        end
    end
    I = setdiff(1:size(mesh.p,2),D)

    # Solution vector
    U = ones(size(mesh.p,2),1)

    ffun = (x,y)->sin(2*pi*x-0.7).*sin(2*pi*y-0.3)

    U=ffun(mesh.p[1,:]',mesh.p[2,:]')

    # Gradient calc
    DX = @bilin_asm(ux*v,mesh)
    DY = @bilin_asm(uy*v,mesh)
    M = @bilin_asm(u*v,mesh)
    CM = cholfact(M)
    dx = w->CM\(DX*w)
    dy = w->CM\(DY*w)

    for i=1:3
        trisurf(mesh,U[:,1])
        X = 1.0/sqrt(1.0+dx(U).^2+dy(U).^2)
        S = @bilin_asm_mat(a*(ux*vx+uy*vy),mesh,X)
        #Solve the linear system
        U[I,1]=cholfact(S[I,I])\(-S[I,D]*U[D,1])
    end

    trisurf(mesh,U[:,1])
end

end
