// C functions used in 
// do not consider censoring group now
// the correpsponding R wrap function is 'crrs.r'

// author: Bingqing Zhou
// Date:   04/22/2010

#include <R.h>
#include <math.h>

// obtain the variance estimate for highly stratified data

void crrvvh(double *t2, int *ici, int *nin, double *x, int *ncov,double *x2, int *ncov2,
      double *ft, int *ndfin, int *strata, int *nstrata, double *wt, double *b, double *v, double *v2)
{
    const int n1=ncov[0], n2=ncov2[0], n=nin[0], ndf=ndfin[0], np=n1+n2, ns=nstrata[0];
	int i,j, k, j1, j2, s, count=0, count2=0, pi[n];
    double a[n][n1],aa[n][n2],tt[ndf][n2], z[np],zb, wye, wyez[np], ss0[n], ss1[n][np], ss2[n][np][np],
		eta[n][np], vt[np][np], q[n][np], xi[ns][np];

	//relates the matrices to the vector parameters
    if (n1 > 0)
        for (i = 0; i < n; i++)
            for (j = 0; j < n1; j++)
                a[i][j]= x[i + n * j];

    if (n2 > 0)
        for (i = 0; i < n; i++)
            for (j = 0; j < n2; j++)
                aa[i][j]= x2[i + n * j];
    if (n2 > 0)
        for (i = 0; i < ndf; i++)
            for (j = 0; j < n2; j++)
                tt[i][j]= ft[i + ndf * j];

    // initialization
	for (i = 0; i < n; i++)
	{
		pi[i] = 0;
		for (j = 0; j < np; j++)
			eta[i][j] = q[i][j] = 0;
	}

	for (i = 0; i < n; i++)
    {
        ss0[i] = 0;
        for (j = 0; j < np; j++)
        {
            ss1[i][j]= 0;
            for (k = 0; k < np; k++)
                ss2[i][j][k] = 0;
        }
    }

	for (i = 0; i < np; i++)
        for (j = 0; j < np; j++)
			vt[i][j] = 0;

	for (i = 0; i < ns; i++)
        for (j = 0; j < np; j++)
			xi[i][j] = 0;

    // obtain s(0)(beta,t),s(1)(beta,t), and s(2)(beta,t) at all type 1 failure times
    for (i = 0; i < n; i++)
    {
        if (ici[i] != 1) continue;
		s = strata[i];

        for (j = 0; j < n; j++)
        {
            if (strata[j] != s) continue;
			if (t2[j] < t2[i] && ici[j] <= 1) continue;
            //put covariates for jth obsertation at Ti in a vector z, and beta'z=zb
            zb = 0;
            for (k = 0; k < np; k++)
                z[k] = 0;

            if (n1 > 0)
                for (j1 = 0; j1 < n1; j1 ++)
                {
                     z[j1] = a[j][j1];
                     zb += b[j1] * a[j][j1];
                }

            if (n2 > 0)
                for (k = 0; k < n2; k++)
                {
                     z[n1+k] = aa[j][k] * tt[count][k];
                     zb += b[n1+k] * aa[j][k] * tt[count][k];
                }

            if (t2[j] >= t2[i])
                wye = exp(zb);
            else
                wye = exp(zb) * wt[i] / wt[j];

            ss0[i] += wye;
            for (k = 0; k < np; k++)
            {
                ss1[i][k] += wye * z[k];
                for (j1 = 0; j1 < np; j1++)
                    ss2[i][k][j1] += wye * z[k] * z[j1];
            }
        }
        count++;
	}

  ///////////////////////////////////////////////////////////////
  
    // eta_i and infomation
	count = 0;
    for (i = 0; i < n; i++)
    {
        for (k = 0; k < np; k++)
            eta[i][k] = 0;

        //dN_i portion of eta_i
        if (ici[i] == 1) 
		{
			zb = 0;
			for (k = 0; k < np; k++)
				z[k] = 0;
			if (n1 > 0)
				for (j = 0; j < n1; j++)
				{
					z[j] = a[i][j];
					zb += b[j] * a[i][j];
				}
			if (n2 > 0)
				for (k = 0; k < n2; k++)
				{
					z[n1+k]= aa[i][k] * tt[count][k];
					zb += b[n1+k] * aa[i][k] * tt[count][k];
				}
      
			for (k = 0; k < np; k++)
				eta[i][k] = z[k] - ss1[i][k]/ss0[i];

			//information
			for (j = 0; j < np; j++)
				for (k = 0; k < np; k++)
					vt[j][k] += ss2[i][j][k]/ss0[i] - ss1[i][j] * ss1[i][k]/ss0[i]/ss0[i];

			count++;
		}
      ////////////////////////////
        // psi, q, pi
        if (ici[i] == 0) 
		{
			pi[i]=0;
			for (j = 0; j < n; j++)
				if (t2[j] >= t2[i])
					pi[i]++;

			count2 = 0;
			for (j1 = 0; j1 < n; j1++)
			{
				if (ici[j1] == 1)  count2 ++;
				if (t2[j1] < t2[i] || ici[j1] != 1) continue;

				//q
				wye = 0;
				for (k = 0; k < np; k++)
					wyez[k] = 0;
				s = strata[j1];
				
				for (j2 = 0; j2 < n; j2++)
				{
					if (t2[j2] >= t2[i]) break;
					if (ici[j2] <= 1 || strata[j2] != s) continue;

					zb = 0;
					if (n1 > 0)
						for (j = 0; j < n1; j++)
						{
							z[j] = a[j2][j];
							zb += b[j] * a[j2][j];
						}
					if (n2 > 0)
						for (k = 0; k < n2; k++)
						{
							z[n1+k]= aa[j2][k] * tt[count2-1][k];
							zb += b[n1+k] * aa[j2][k] * tt[count2-1][k];
						}

					wye += exp(zb)* wt[j1] / wt[j2];
					for (k = 0; k < np; k++)
						wyez[k] += exp(zb)* wt[j1] / wt[j2] * z[k];
				}

				for (k = 0; k < np; k++)
					q[i][k] += (wyez[k] - ss1[j1][k]/ss0[j1] * wye) /ss0[j1];
			}
		}
	}

  ///////////////////////////////////////////////////////////////
    //result
	for (i = 0; i < n; i++)
        for (k = 0; k < np; k++)
		{
			if (ici[i]==0) eta[i][k] += q[i][k]/pi[i];
		    for (j = 0; j < n; j++)
			{
				    if (t2[j] > t2[i]) break;
				    if (ici[j] == 0) eta[i][k] -= q[j][k]/pi[j]/pi[j];
			 }
		}

	//output information
    for (j = 0; j < np; j++)
        for (k= 0; k < np; k++)
            v[j  + k * np] = vt[j][k] ;   

	//output variance of the estimating equation
	for (i = 0; i < np*np; i++)
		v2[i] = 0;
	for (i = 0; i < ns; i++)
		for (j = 0; j < n; j++)
		{
			if (strata[j] != i+1) continue;
			for (k = 0; k < np; k++)
				xi[i][k] += eta[j][k];
		}
	for (j = 0; j < np; j++)
		for (k= 0; k < np; k++)
			for (i = 0; i < ns; i++)
				v2[j  + k * np] += xi[i][j] * xi[i][k];
}




// reproduce crrf to obtain liklihood

void crrf(double *t2, int *ici, int *nin, double *x, int *ncov,  double *x2, int *ncov2, 
		 double *ft, int *ndfin, double *wt, double *b, double *lik)
{
    int i,j, k, j1,  count=0;
	const int n1=ncov[0], n2=ncov2[0], n=nin[0], ndf=ndfin[0];
	double a[n][n1], aa[n][n2], tt[ndf][n2], likli=0, ss0,zb;
	
	if (n1 > 0)
        for (i = 0; i < n; i++)
            for (j = 0; j < n1; j++)
                a[i][j]= x[i + n * j];

    if (n2 > 0)
        for (i = 0; i < n; i++)
            for (j = 0; j < n2; j++)
                aa[i][j]= x2[i + n * j];
    if (n2 > 0)
        for (i = 0; i < ndf; i++)
            for (j = 0; j < n2; j++)
                tt[i][j]= ft[i + ndf * j];

    for (i = 0; i < n; i++)
	{
		if (ici[i] != 1) continue;

		//first part of logliklihood
		if (n1 > 0)
			for (j = 0; j < n1; j++)
				likli += b[j] * a[i][j];
		if (n2 > 0)
			for (k = 0; k < n2; k++)
				likli += b[n1+k] * aa[i][k] * tt[count][k];
        //second part of logliklihood
		ss0=0;
		for (j = 0; j < n; j++)
		{
			if (t2[j] < t2[i] && ici[j] <= 1) continue;
			zb = 0.0;
			if (n1 > 0)
				for (j1 = 0; j1 < n1; j1 ++)
					zb += b[j1] * a[j][j1];
			if (n2 > 0)
				for (k = 0; k < n2; k++)
					zb += b[n1+k] * aa[j][k] * tt[count][k];

			if (t2[j] >= t2[i]) 
				ss0 += exp(zb);
			else 
				ss0 += exp(zb) * wt[i] / wt[j]; 

		}

		likli -= log(ss0);

	    count++;
	}

   
	*lik=likli;

}


// reproduce crrfsv to obtain liklihood, score, and information

void crrfsv(double *t2, int *ici, int *nin, double *x, int *ncov,  double *x2, int *ncov2, 
		 double *ft, int *ndfin, double *wt, double *b, double *lik, double *st, double *v)
{
    int i,j, k, j1,  count=0;
	const int n1=ncov[0], n2=ncov2[0], n=nin[0], ndf=ndfin[0], np=n1+n2;
	double a[n][n1],aa[n][n2],tt[ndf][n2],likli=0,ss0,ss1[np],ss2[np][np],z[np],zb,wye,vt[np][np];

    for (i = 0; i < np; i++)
	{
		st[i] = 0;
		for (j = 0; j < np; j++)
			vt[i][j]= 0;
	}

	if (n1 > 0)
		for (i = 0; i < n; i++)
			for (j = 0; j < n1; j++)
				a[i][j]= x[i + n * j];

    if (n2 > 0)
		for (i = 0; i < n; i++)
		    for (j = 0; j < n2; j++)
			    aa[i][j]= x2[i + n * j];
    if (n2 > 0)
		for (i = 0; i < ndf; i++)
		    for (j = 0; j < n2; j++)   
		     	tt[i][j]= ft[i + ndf * j];


    for (i = 0; i < n; i++)
	{
		if (ici[i] != 1) continue;

		//first part of logliklihood and score
		zb = 0;
		if (n1 > 0)
			for (j = 0; j < n1; j++)
			{
				st[j]  +=  a[i][j];
				zb += b[j] * a[i][j];
			}
		if (n2 > 0)
			for (k = 0; k < n2; k++)
			{
				st[n1+k]  +=  aa[i][k] * tt[count][k];
				zb += b[n1+k] * aa[i][k] * tt[count][k];
			}

        likli += zb;

        //second part of logliklihood, second part of score, and information
		ss0 = 0;
		for (k = 0; k < np; k++)
		{
			ss1[k] = 0;
			for (j1 = 0; j1 < np; j1 ++)
				ss2[k][j1] = 0;
		}

		for (j = 0; j < n; j++)
		{
			if (t2[j] < t2[i] && ici[j] <= 1) continue;

			zb = 0;
			for (k = 0; k < np; k++)
				z[k] = 0;

			if (n1 > 0)
				for (j1 = 0; j1 < n1; j1 ++)
				{
					z[j1] = a[j][j1];
					zb += b[j1] * a[j][j1];
				}
				    
			if (n2 > 0)
				for (k = 0; k < n2; k++)
				{
					z[n1+k] = aa[j][k] * tt[count][k];
					zb += b[n1+k] * aa[j][k] * tt[count][k];
				}


			if (t2[j] >= t2[i]) 
				wye = exp(zb);
			else 
				wye = exp(zb) * wt[i] / wt[j]; 

			ss0 += wye;
			for (k = 0; k < np; k++)
			{
				ss1[k] += wye * z[k];
				for (j1 = 0; j1 < np; j1++)
					ss2[k][j1] += wye * z[k] * z[j1];
			}

		}

		likli -= log(ss0);

		for (k = 0; k < np; k++)
		{
			st[k] -= ss1[k] / ss0;
			for (j1 = 0; j1 < np; j1++)
				vt[k][j1] += ss2[k][j1]/ss0 - ss1[k] * ss1[j1]/ss0/ss0;
		}
		
	    count++;
	}

	*lik=likli;

    for (i = 0; i < np; i++)
		for (j = 0; j < np; j++)
			v[i + np * j] = vt[i][j];

}





  // reproduce crrvv to obtain the variance 

void crrvv(double *t2, int *ici, int *nin, double *x, int *ncov,double *x2, int *ncov2,
         double *ft, int *ndfin, double *wt, double *b, double *v, double *v2)
{
    const int n1=ncov[0], n2=ncov2[0], n=nin[0], ndf=ndfin[0], np=n1+n2;
	int i,j, k, j1, j2,count=0, count1=0, count2=0, pi[n];
    double a[n][n1],aa[n][n2],tt[ndf][n2], z[np],zb, wye, wyez[np], ss0[n], ss1[n][np], ss2[n][np][np],
		eta[n][np], vt[np][np], q[n][np];

    // initialization
	for (i = 0; i < n; i++)
	{
		pi[i] = 0;
		for (j = 0; j < np; j++)
			eta[i][j] = q[i][j] = 0;
	}

	for (i = 0; i < n; i++)
    {
        ss0[i] = 0;
        for (j = 0; j < np; j++)
        {
            ss1[i][j]= 0;
            for (k = 0; k < np; k++)
                ss2[i][j][k] = 0;
        }
    }

	for (i = 0; i < np; i++)
        for (j = 0; j < np; j++)
			vt[i][j] = 0;

    if (n1 > 0)
        for (i = 0; i < n; i++)
            for (j = 0; j < n1; j++)
                a[i][j]= x[i + n * j];

    if (n2 > 0)
        for (i = 0; i < n; i++)
            for (j = 0; j < n2; j++)
                aa[i][j]= x2[i + n * j];
    if (n2 > 0)
        for (i = 0; i < ndf; i++)
            for (j = 0; j < n2; j++)
                tt[i][j]= ft[i + ndf * j];

    // obtain s(0)(beta,t),s(1)(beta,t), and s(2)(beta,t) at all type 1 failure times
    for (i = 0; i < n; i++)
    {
        if (ici[i] != 1) continue;

        for (j = 0; j < n; j++)
        {
            if (t2[j] < t2[i] && ici[j] <= 1) continue;
            //put covariates for jth obsertation at Ti in a vector z, and beta'z=zb
            zb = 0;
            for (k = 0; k < np; k++)
                z[k] = 0;

            if (n1 > 0)
                for (j1 = 0; j1 < n1; j1 ++)
                {
                    z[j1] = a[j][j1];
                    zb += b[j1] * a[j][j1];
                }

            if (n2 > 0)
                for (k = 0; k < n2; k++)
                {
                    z[n1+k] = aa[j][k] * tt[count][k];
                    zb += b[n1+k] * aa[j][k] * tt[count][k];
                }

            if (t2[j] >= t2[i])
                wye = exp(zb);
            else
                wye = exp(zb) * wt[i] / wt[j];

            ss0[i] += wye;
            for (k = 0; k < np; k++)
            {
                ss1[i][k] += wye * z[k];
                for (j1 = 0; j1 < np; j1++)
                    ss2[i][k][j1] += wye * z[k] * z[j1];
            }
        }
        count++;
    }

  ///////////////////////////////////////////////////////////////
  
    count = 0;
    // eta_i and infomation
    for (i = 0; i < n; i++)
    {
        for (k = 0; k < np; k++)
            eta[i][k] = 0;

		count1 = 0;
        //dLamda portion of eta_i
        for (j = 0; j < n; j++)
        {
            if (ici[j] != 1) continue;

            zb = 0;
            for (k = 0; k < np; k++)
                z[k] = 0;

            if (n1 > 0)
                for (j1 = 0; j1 < n1; j1 ++)
                {
                    z[j1] = a[i][j1];
                    zb += b[j1] * a[i][j1];
                }

            if (n2 > 0)
                for (k = 0; k < n2; k++)
                {
                    z[n1+k] = aa[i][k] * tt[count1][k];
                    zb += b[n1+k] * aa[i][k] * tt[count1][k];
                }


            if (t2[j] <= t2[i])
                wye = exp(zb);
            else if (t2[j] > t2[i] && ici[i]>1)
                wye = exp(zb) * wt[j] / wt[i];
            else continue;

            for (k = 0; k < np; k++)
                eta[i][k] -= (z[k] - ss1[j][k]/ss0[j]) * wye / ss0[j];

            count1 ++;
        }

        //dN_i portion of eta_i
        if (ici[i] == 1) 
		{
			zb = 0;
			for (k = 0; k < np; k++)
				z[k] = 0;
			if (n1 > 0)
				for (j = 0; j < n1; j++)
				{
					z[j] = a[i][j];
					zb += b[j] * a[i][j];
				}
			if (n2 > 0)
				for (k = 0; k < n2; k++)
				{
					z[n1+k]= aa[i][k] * tt[count][k];
					zb += b[n1+k] * aa[i][k] * tt[count][k];
				}

			for (k = 0; k < np; k++)
				eta[i][k] += z[k] - ss1[i][k]/ss0[i];

			//information
			for (j = 0; j < np; j++)
				for (k = 0; k < np; k++)
					vt[j][k] += ss2[i][j][k]/ss0[i] - ss1[i][j] * ss1[i][k]/ss0[i]/ss0[i];

			count++;
		}
      ////////////////////////////
        // psi, q, pi
        if (ici[i] == 0) 
		{
			pi[i]=0;
			for (j = 0; j < n; j++)
				if (t2[j] >= t2[i])
					pi[i]++;

			count2 = 0;
			for (j1 = 0; j1 < n; j1++)
			{
				if (ici[j1] == 1)  count2 ++;
				if (t2[j1] < t2[i] || ici[j1] != 1) continue;

				//q
				wye = 0;
				for (k = 0; k < np; k++)
					wyez[k] = 0;
				
				for (j2 = 0; j2 < n; j2++)
				{
					if (t2[j2] >= t2[i]) break;
					if (ici[j2] <= 1) continue;


					zb = 0;
					if (n1 > 0)
						for (j = 0; j < n1; j++)
						{
							z[j] = a[j2][j];
							zb += b[j] * a[j2][j];
						}
					if (n2 > 0)
						for (k = 0; k < n2; k++)
						{
							z[n1+k]= aa[j2][k] * tt[count2-1][k];
							zb += b[n1+k] * aa[j2][k] * tt[count2-1][k];
						}

					wye += exp(zb)* wt[j1] / wt[j2];
					for (k = 0; k < np; k++)
						wyez[k] += exp(zb)* wt[j1] / wt[j2] * z[k];
				}

				for (k = 0; k < np; k++)
					q[i][k] += (wyez[k] - ss1[j1][k]/ss0[j1] * wye) /ss0[j1];
			}
		}
	}

  ///////////////////////////////////////////////////////////////
    //output result
	for (i = 0; i < n; i++)
        for (k = 0; k < np; k++)
		{
			if (ici[i]==0) eta[i][k] += q[i][k]/pi[i];
		    for (j = 0; j < n; j++)
			{
				    if (t2[j] > t2[i]) break;
				    if (ici[j] == 0) eta[i][k] -= q[j][k]/pi[j]/pi[j];
			 }
		}

    for (j = 0; j < np; j++)
        for (k= 0; k < np; k++)
            v[j  + k * np] = vt[j][k] ;   

	for (i=0; i<np*np; i++)
		v2[i]=0;
    
    for (j = 0; j < np; j++)
        for (k= 0; k < np; k++)
			for (i = 0; i < n; i++)
                v2[j  + k * np] += eta[i][j] * eta[i][k] ;

}


//end of the functions

// example
/*

setwd( "C:/fine/practice")
dyn.load("crrc.dll")
library(cmprsk)
set.seed(1671)
n=200
np1=3
ftime <- sort(rexp(n))
fstatus <- sample(0:2,n,replace=TRUE)
cov <- matrix(runif(np1*n),nrow=n)
cov2 = cbind(cov[,1],cov[,1])
np = np1 + ncol(cov2)
npt = ncol(cov2)
uft <- ftime[fstatus==1]
tfs=cbind(uft,uft^2)
ndf=length(uft)
b= c(-.2,.7,1,-.5,3)
cengroup = rep(1,length(fstatus))
ncg = 1
cenind <- ifelse(fstatus==1,1,0)
uuu <- matrix(0,nrow=ncg,ncol=length(ftime))

for (k in 1:ncg) {
    u <- do.call('survfit',list(formula=Surv(ftime,cenind)~1,data=data.frame(ftime,cenind,cengroup),subset=cengroup==k))
    ### note: want censring dist km at ftime-
    u <- summary(u,times=sort(ftime*(1-.Machine$double.eps)))
    uuu[k,1:length(u$surv)] <- u$surv
}

  #####################
  start1 <- Sys.time()
.C("crrf", as.double(ftime),as.integer(fstatus),as.integer(length(ftime)),
   as.double(cov), as.integer(np-npt), as.double(cov2),as.integer(npt),
   as.double(tfs),as.integer(ndf), as.double(uuu), 
   as.double(b), double(1))[[12]]
start2 <- Sys.time()
.Fortran('crrf',as.double(ftime),as.integer(fstatus),
                  as.integer(length(ftime)),as.double(cov),as.integer(np-npt),
                  as.integer(np),as.double(cov2),as.integer(npt),
                  as.double(tfs),as.integer(ndf),as.double(uuu),
                  as.integer(ncg),as.integer(cengroup),as.double(b),
                  double(1),double(np),package="cmprsk")[15:16]
end2 <- Sys.time()
start2-start1
end2-start2

  ####################

    start1 = Sys.time()
  .C("crrfsv", as.double(ftime),as.integer(fstatus),as.integer(length(ftime)),
   as.double(cov), as.integer(np-npt), as.double(cov2),as.integer(npt),
   as.double(tfs),as.integer(ndf), as.double(uuu), 
   as.double(b), double(1), double(np), double(np*np))[12:14]

  start2 = Sys.time()
  .Fortran('crrfsv',as.double(ftime),as.integer(fstatus),
                  as.integer(length(ftime)),as.double(cov),as.integer(np-npt),
                  as.integer(np),as.double(cov2),as.integer(npt),
                  as.double(tfs),as.integer(ndf),as.double(uuu),
                  as.integer(ncg),as.integer(cengroup),as.double(b),
                  double(1),double(np),double(np*np),double(np),double(np),
                  double(np*np))[15:17]
   end = Sys.time()
    
   start2-start1
   end-start2
  
	#############################


  start <- Sys.time()
.C("crrvv", as.double(ftime),as.integer(fstatus),as.integer(length(ftime)),
   as.double(cov), as.integer(np-npt), as.double(cov2),as.integer(npt),as.double(tfs),
   as.integer(ndf), as.double(uuu),  as.double(b), double(np*np), double(np*np))[12:13]
  start1 <- Sys.time()
  
 # uft <- unique(ftime[fstatus==1])
 # tfs=cbind(uft,uft^2)
 # ndf=length(uft)
   .Fortran('crrvv',as.double(ftime),as.integer(fstatus),
                as.integer(length(ftime)),as.double(cov),as.integer(np-npt),
                as.integer(np),as.double(cov2),as.integer(npt),
                as.double(tfs),as.integer(ndf),as.double(uuu),
                as.integer(ncg),as.integer(cengroup),as.double(b),
                double(np*np),double(np*np),double(np*np),
                double(length(ftime)*(np+1)),double(np),double(np),
                double(2*np),double(np),double(ncg*np),integer(ncg), package="cmprsk")[15:16]
   end <- Sys.time()
   start1-start
   end-start1
*/



