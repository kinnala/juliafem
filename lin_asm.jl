macro lin_asm_dispatch(lin_form,mesh,nodal...)
    quote
        ##
        ## Assembly of a linear form using P1 elements.
        ## Results in a vector f_i.
        ##
        local I = Array(Int64,0)
        local V = Array(Float64,0)

        # Local basis functions and their gradients
        local phi_hat = [(x,y)->1-x-y, (x,y)->x, (x,y)->y]
        local dphi_hat = [-1.0 -1.0; 1.0 0.0; 0.0 1.0]'

        # Dunavant quadrature rule. This is in fact
        # good for 2nd order polynomials on triangle.
        local qp = [0.5 0; 0.5 0.5; 0 0.5]'
        local qw = [1/6,1/6,1/6]
        local Nqp = size(qp,2)

        # Array for storing the pre-evaluated basis functions
        # at the quadrature points
        local phi_hat_qp = Array(Float64,3,Nqp)

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
        local A = Array(Float64,2,2)
        local invAt = Array(Float64,2,2)
        local n1 = Array(Float64,2)
        local n2 = Array(Float64,2)
        local n3 = Array(Float64,2)
        local n1ix::Int64 = -1
        local n2ix::Int64 = -1
        local n3ix::Int64 = -1
        local gqp = Array(Float64,2,Nqp)
        local dphi = Array(Float64,size(dphi_hat,1),size(dphi_hat,2))
        local absdetA::Float64 = -1.0
        local detA::Float64 = -1.0
        if length($nodal)==1
            local a1::Float64 = -1.0
            local a2::Float64 = -1.0
            local a3::Float64 = -1.0
            local aset = $(length(nodal)==1 ? nodal[1] : nothing)
        end

        for k=1:size($mesh.t,2)
            # Find global node indices
            n1ix = $mesh.t[1,k]
            n2ix = $mesh.t[2,k]
            n3ix = $mesh.t[3,k]
            # Find the parameter values at nodes
            if length($nodal)==1
                a1 = aset[n1ix,1]
                a2 = aset[n2ix,1]
                a3 = aset[n3ix,1]
            end
            # Find nodes
            n1[1] = $mesh.p[1,n1ix]
            n1[2] = $mesh.p[2,n1ix]
            n2[1] = $mesh.p[1,n2ix]
            n2[2] = $mesh.p[2,n2ix]
            n3[1] = $mesh.p[1,n3ix]
            n3[2] = $mesh.p[2,n3ix]
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
                local vx::Float64 = dphi[1,i]
                local vy::Float64 = dphi[2,i]
                local ix::Int64 = $mesh.t[i,k]
                # Loop over quadrature points
                for q=1:Nqp
                    x = gqp[1,q]
                    y = gqp[2,q]
                    v = phi_hat_qp[i,q]
                    # TODO try to calculate these only once
                    if length($nodal)==1
                        a = a1+(a2-a1)*qp[1,q]+(a3-a1)*qp[2,q]
                    end
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
