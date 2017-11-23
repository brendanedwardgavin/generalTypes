#
function bicgstab(A,x0,b,maxit,epsmax)
    x=copy(x0)
    r=b-A*x
    rtilde=copy(r)
    p=copy(r)

    for i in 0:maxit-1
        #solve M*phat=p
        phat=copy(p)
        v=A*phat
        alpha=(r'*rtilde)[1]/(v'*rtilde)[1]
        r2=r-alpha*v
        #println("$i $(norm(r2)/norm(b))")
        if(norm(r2)<epsmax)
            x=x+alpha*phat
            return x
        end

        #solve M*rhat=r2
        rhat=copy(r2)
        u=A*rhat
        omega=(r2'*u)[1]/(u'*u)[1]
        x=x+alpha*phat+omega*rhat
        rold=copy(r)
        r=r2-omega*u
        println("$i $(norm(r)/norm(b))")
        if(norm(r)<epsmax)
            return x
        end

        beta=(alpha/omega)*((r'*rtilde)[1]/(rold'*rtilde)[1])
        p=r+beta*(p-omega*v)
    end

    return x
end



function zbicgstab(A,b,x0,maxit,eps)
	#matrix multiply:

        r=b-A*x0
	s=copy(r)
	y0=copy(r)

	x=copy(x0)

	#inner product:
	delta=(y0'*r)[1]
	#phi=(y0'*(A*s))[1]/delta

	for i in 1:maxit

		#matrix multiply:
		As=A*s
		#inner product
		phi=(y0'*(As))[1]/delta
		omega=1/phi
		w=r-As*omega

		#matrix multiply:
		Aw=A*w
		#inner product:
		chi=((Aw)'*w)[1]/norm(Aw)^2
		r=w-chi*Aw


		#println("   $i     $(norm(r))   $(norm(Aw))")
		if(norm(r)<eps)
			#println("r finish: $(norm(r))")
			return x
		end

		x=x+s*omega+w*chi
		deltaold=delta
		#inner product:
		delta=(y0'*r)[1]
		psi=-1.0*omega*delta/(deltaold*chi)
		s=r-(s-(As)*chi)*psi
		#print("     $(norm(s))\n")
		if(norm(s)<eps)
			#println("s finish: $(norm(s))")
			return x
		end


	end

	#println("          no converge: $(norm(r))")
	return x
end
