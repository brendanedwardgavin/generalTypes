


function trapezoidal(nc)
    points=zeros(nc)
    weights=zeros(nc)

    a=-1.0
    b=1.0
    dx=(b-a)/(nc)
    for i in 1:nc
        points[i]=a+dx*(i-1)
        weights[i]=(b-a)/(nc+1)
    end

    return (points,weights)
end



#linear feast for finding right eigenvectors
function feast_nsR(A,B,x0,nc,emid,r,eps,maxit)

    (n,m0)=size(x0)
    lest=zeros(m0)
    x=copy(x0)
   
    (gk,wk)=trapezoidal(nc)
    offset=pi/nc #offset angle of first quadrature point; makes sure it isn't on real axis

    residuals=zeros(m0)
    p=sortperm(residuals)

    it=0
    res=1.0
    ninside=m0
    res=1.0
    Q=copy(x)
    Qk=zeros(x)
    while it<maxit && abs.(res)>abs.(eps)
        it=it+1 

        #orthogonalize subspace with SVD:
        Bq=Q'*B*Q
        (ll,qq)=eig(Bq)
        
        #println(abs.(ll)/maximum(abs.(ll)))
        nremove=0
        maxll=maximum(abs.(ll))
        for i in 1:m0
            if(abs.(ll[i]/maxll)<1e-10) #arbitrary singular value cutoff
                nremove=nremove+1
            end
        end
        if(nremove>0)
            println("Subspace too large, resizing:",m0-nremove)
        end
        println("     Smallest s=",minimum(abs.(ll[nremove+1:m0]))/maximum(abs.(ll)),)
        p=sortperm(abs.(ll))
        U=Q*qq[:,p[nremove+1:m0]]*diagm(1./sqrt.(ll[p[nremove+1:m0]]))
        m0=m0-nremove
        #(Qq,Rq)=qr(Q)
        
        Aq=U'*A*U
        Bq=U'*B*U

        (lest,xq)=eig(Aq,Bq)

        x=U*xq

        for j in 1:m0
            x[:,j]=x[:,j]/norm(x[:,j])
        end
        
        resvecs=A*x-B*x*diagm(lest)
        residuals=zeros(m0)

        for j in 1:m0
            residuals[j]=norm(resvecs[:,j])
        end 

        ninside=0
        ninside_real=0
        for j in 1:m0
            if(abs.(lest[j]-emid)<abs.(r))
                ninside=ninside+1
            end
        end

        if(ninside>m0)
            error("Number of eigenvalues inside contour is greater than the subspace dimension M0")
        end

        ninside_real=ninside
        if(ninside==0)
            ninside = m0
        end
        
        #println(lest," ",ninside_real)

        #get distances of estimated eigenvalues from middle of contour:
        lestdif=abs.(lest.-emid)
        #sort eigenpairs based on those distances:
        p=sortperm(lestdif)     
 
        #p=sortperm(residuals)
        res=residuals[p[ninside]]
        
        println("    $it: $(residuals[p[ninside]])   $ninside_real")

        if(res<eps)
            break
        end

        #=Q[:]=0.0
        for j in 1:nc
            z=emid+r*exp(im*(gk[j]*pi+pi+offset))
            #Qk=\(z*B-A,B*x)
            Qk[:]=0.0
            for i in 1:m0
                Qk[:,i]=\(z*B-A,resvecs[:,i]) 
                #Qk[:,i]=zbicgstab(z*B-A,resvecs[:,i],zeros(n,1),100,1e-2)
            end
            Q=Q+wk[j]*exp(im*(gk[j]*pi+pi+offset))*(x-Qk)*diagm(1.0./(z.-lest))
        end=#

        Q[:]=0.0
        #bcgx0=zeros(Q[:,1])
	bcgx0=zeros(Q)
        for j in 1:nc
            println("        CP $j of $nc")
            z=emid+r*exp(im*(gk[j]*pi+pi+offset))
            #Qk=\(z*B-A,B*x)
            #Qk[:]=0.0
            Bx=B*x
            Qk=zbicgstabBlock(z*B-A,Bx,bcgx0,500,min(res*1e-2,1e-2))
            #for i in 1:m0
                #Qk[:,i]=\(z*B-A,(B*x)[:,i])
            #    Qk[:,i]=zbicgstab(z*B-A,Bx[:,i],bcgx0,500,min(res*1e-2,1e-2))
            #end
            Q=Q+wk[j]*exp(im*(gk[j]*pi+pi+offset))*Qk
        end
    end
   
   return (lest,x) 
end



#linear feast for finding right eigenvectors
function pfeast_nsR(A,B,x0,nc,emid,r,eps,maxit,procs)

    sprocs=sort(procs)

    (n,m0)=size(x0)
    lest=zeros(m0)
    x=copy(x0)
   
    (gk,wk)=trapezoidal(nc)
    offset=pi/nc #offset angle of first quadrature point; makes sure it isn't on real axis

    residuals=zeros(m0)
    p=sortperm(residuals)

    it=0
    res=1.0
    ninside=m0
    res=1.0
    Q=copy(x)
    Qk=zeros(x)
    while it<maxit && abs.(res)>abs.(eps)
        it=it+1 

        #orthogonalize subspace with SVD:
        Bq=Q'*B*Q
        (ll,qq)=eig(Bq)
        U=Q*qq*diagm(1./sqrt.(ll))
        #(Qq,Rq)=qr(Q)
        println("     Smallest s=",minimum(abs.(ll)))
        Aq=U'*A*U
        Bq=U'*B*U

        (lest,xq)=eig(Aq,Bq)

        x=U*xq

        for j in 1:m0
            x[:,j]=x[:,j]/norm(x[:,j])
        end
        
        resvecs=A*x-B*x*diagm(lest)
        residuals=zeros(m0)

        for j in 1:m0
            residuals[j]=norm(resvecs[:,j])
        end 

        ninside=0
        ninside_real=0
        for j in 1:m0
            if(abs.(lest[j]-emid)<abs.(r))
                ninside=ninside+1
            end
        end

        if(ninside>m0)
            error("Number of eigenvalues inside contour is greater than the subspace dimension M0")
        end

        ninside_real=ninside
        if(ninside==0)
            ninside = m0
        end
        
        #println(lest," ",ninside_real)

        #get distances of estimated eigenvalues from middle of contour:
        lestdif=abs.(lest.-emid)
        #sort eigenpairs based on those distances:
        p=sortperm(lestdif)     
 
        #p=sortperm(residuals)
        res=residuals[p[ninside]]
        
        println("    $it: $(residuals[p[ninside]])   $ninside_real")

        if(res<eps)
            break
        end

        #=Q[:]=0.0
        for j in 1:nc
            z=emid+r*exp(im*(gk[j]*pi+pi+offset))
            #Qk=\(z*B-A,B*x)
            Qk[:]=0.0
            for i in 1:m0
                Qk[:,i]=\(z*B-A,resvecs[:,i]) 
                #Qk[:,i]=zbicgstab(z*B-A,resvecs[:,i],zeros(n,1),100,1e-2)
            end
            Q=Q+wk[j]*exp(im*(gk[j]*pi+pi+offset))*(x-Qk)*diagm(1.0./(z.-lest))
        end=#

        #Q[:]=0.0
        bcgx0=zeros(Q[:,1])
        results=Vector{Future}(nc-1) 

        for j in 2:nc
            results[j-1]=@spawnat sprocs[j] solveFeastCP(j,nc,emid,r,offset,res,wk,gk,A,B,x) 
        end
        j=1
        Q=solveFeastCP(j,nc,emid,r,offset,res,wk,gk,A,B,x)
        
        for j in 2:nc
            wait(results[j-1])
        end

        for j in 2:nc
            Q=Q+fetch(results[j-1])
        end
    end
   
   return (lest,x) 
end

function solveFeastCP(k,nc,emid,r,offset,res,wk,gk,A,B,x)
    m0=size(x,2)
    Qk=zeros(x)
    println("        CP $k of $nc")
    z=emid+r*exp(im*(gk[k]*pi+pi+offset))
    #Qk=\(z*B-A,B*x)
    Bx=B*x
    bcgx0=zeros(Qk[:,1])
    for i in 1:m0
        #Qk[:,i]=\(z*B-A,(B*x)[:,i])
        Qk[:,i]=zbicgstab(z*B-A,Bx[:,i],bcgx0,500,min(res*1e-2,1e-2))
    end
    return (wk[k]*exp(im*(gk[k]*pi+pi+offset)))*Qk
end
