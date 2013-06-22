// C functions for clustered fine-gray model
// for simplicity, do not consider censoring group now
// the correpsponding R wrap function is 'crrcc.r'

// author: Bingqing Zhou
// Date:   12/21/2008
// modified on 5/19/2013

#include <R.h>
#include <math.h>
#include <stdlib.h>
#define LEN sizeof(double)

// similar to crrf to obtain liklihood

void crrfoo(double *t2, int *ici, int *nin, double *x, int *ncov,  double *x2, int *ncov2, 
		 double *ft, int *ndfin, double *wt, double *b, double *lik)
{
    int i,j, k, j1,  count=0;
	const int n1=ncov[0], n2=ncov2[0], n=nin[0], ndf=ndfin[0];
	double likli=0, s0,zb,*a,*aa,*tt;

	a=(double*)malloc(n*n1*LEN);
	aa=(double*)malloc(n*n2*LEN);
	tt=(double*)malloc(ndf*n2*LEN);
	
	if (n1 > 0)
        for (i = 0; i < n; i++)
            for (j = 0; j < n1; j++)
                *(a+i*n1+j)= x[i + n * j];

    if (n2 > 0)
        for (i = 0; i < n; i++)
            for (j = 0; j < n2; j++)
                *(aa+i*n2+j)= x2[i + n * j];
    if (n2 > 0)
        for (i = 0; i < ndf; i++)
            for (j = 0; j < n2; j++)
                *(tt+i*n2+j)= ft[i + ndf * j];

    for (i = 0; i < n; i++)
	{
		if (ici[i] != 1) continue;

		//first part of logliklihood
		if (n1 > 0)
			for (j = 0; j < n1; j++)
				likli += b[j] *  (*(a+i*n1+j)) ;
		if (n2 > 0)
			for (k = 0; k < n2; k++)
				likli += b[n1+k] * (*(aa+i*n2+k))* (*(tt+count*n2+k));
        //second part of logliklihood
		s0=0;
		for (j = 0; j < n; j++)
		{
			if (t2[j] < t2[i] && ici[j] <= 1) continue;
			zb = 0.0;
			if (n1 > 0)
				for (j1 = 0; j1 < n1; j1 ++)
					zb += b[j1] * (*(a+j*n1+j1));
			if (n2 > 0)
				for (k = 0; k < n2; k++)
					zb += b[n1+k] * (*(aa+j*n2+k)) * (*(tt+count*n2+k));

			if (t2[j] >= t2[i]) 
				s0 += exp(zb);
			else 
				s0 += exp(zb) * wt[i] / wt[j]; 

		}

		likli -= log(s0);

	    count++;
	}

   
	*lik=likli;

	free(a);free(aa);free(tt);
}


// reproduce crrfsv to obtain liklihood, score, and information

void crrfsvoo(double *t2, int *ici, int *nin, double *x, int *ncov,  double *x2, int *ncov2, 
		 double *ft, int *ndfin, double *wt, double *b, double *lik, double *st, double *v)
{
    int i,j, k, j1,  count=0;
	const int n1=ncov[0], n2=ncov2[0], n=nin[0], ndf=ndfin[0], np=n1+n2;
	double likli=0,s0,s1[np],z[np],zb,wye,*a,*aa,*tt,*s2,*vt;

	a=(double*)malloc(n*n1*LEN);
	aa=(double*)malloc(n*n2*LEN);
	tt=(double*)malloc(ndf*n2*LEN);
	s2=(double*)malloc(np*np*LEN);
	vt=(double*)malloc(np*np*LEN);

    for (i = 0; i < np; i++)
	{
		st[i] = 0;
		for (j = 0; j < np; j++)
			*(vt+i*np+j)= 0;
	}

	if (n1 > 0)
		for (i = 0; i < n; i++)
			for (j = 0; j < n1; j++)
				*(a+i*n1+j)= x[i + n * j];

    if (n2 > 0)
		for (i = 0; i < n; i++)
		    for (j = 0; j < n2; j++)
			    *(aa+i*n2+j)= x2[i + n * j];
    if (n2 > 0)
		for (i = 0; i < ndf; i++)
		    for (j = 0; j < n2; j++)   
		     	*(tt+i*n2+j)= ft[i + ndf * j];


    for (i = 0; i < n; i++)
	{
		if (ici[i] != 1) continue;

		//first part of logliklihood and score
		zb = 0;
		if (n1 > 0)
			for (j = 0; j < n1; j++)
			{
				st[j]  +=  *(a+i*n1+j);
				zb += b[j] * (*(a+i*n1+j));
			}
		if (n2 > 0)
			for (k = 0; k < n2; k++)
			{
				st[n1+k]  +=  (*(aa+i*n2+k))* (*(tt+count*n2+k));
				zb += b[n1+k] * (*(aa+i*n2+k)) * (*(tt+count*n2+k));
			}

        likli += zb;

        //second part of logliklihood, second part of score, and information
		s0 = 0;
		for (k = 0; k < np; k++)
		{
			s1[k] = 0;
			for (j1 = 0; j1 < np; j1 ++)
				*(s2+k*np+j1)= 0;
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
					z[j1] = *(a+j*n1+j1);
					zb += b[j1] * (*(a+j*n1+j1));
				}
				    
			if (n2 > 0)
				for (k = 0; k < n2; k++)
				{
					z[n1+k] =  (*(aa+j*n2+k)) * (*(tt+count*n2+k));
					zb += b[n1+k] *  (*(aa+j*n2+k)) * (*(tt+count*n2+k));
				}


			if (t2[j] >= t2[i]) 
				wye = exp(zb);
			else 
				wye = exp(zb) * wt[i] / wt[j]; 

			s0 += wye;
			for (k = 0; k < np; k++)
			{
				s1[k] += wye * z[k];
				for (j1 = 0; j1 < np; j1++)
					*(s2+k*np+j1) += wye * z[k] * z[j1];
			}

		}

		likli -= log(s0);

		for (k = 0; k < np; k++)
		{
			st[k] -= s1[k] / s0;
			for (j1 = 0; j1 < np; j1++)
				*(vt+k*np+j1) += (*(s2+k*np+j1))/s0 - s1[k] * s1[j1]/s0/s0;
		}
		
	    count++;
	}

	*lik=likli;

    for (i = 0; i < np; i++)
		for (j = 0; j < np; j++)
			v[i + np * j] = *(vt+i*np+j);


	free(a);free(aa);free(tt);free(s2);free(vt);
}





  // reproduce crrvv to obtain the variance
  // Consider cluster (can be used when # of cluster=1)

void crrvvoo(double *t2, int *ici, int *nin, double *x, int *ncov,double *x2, int *ncov2,
         double *ft, int *nfin, int *cluster, int *ncin, double *wt, double *b, double *v, double *v2)
{
    const int n1=ncov[0], n2=ncov2[0], n=nin[0], ndf=nfin[0], np=n1+n2, nc=*ncin;
	int i,j, k, j1, j2,count=0, count1=0, count2=0, pi[n];
    double z[np],zb, wye, wyez[np], ss0[n],*a,*aa,*tt,*ss1,*ss2,*eta,*q,*vt,*xi;


	a=(double*)malloc(n*n1*LEN);
	aa=(double*)malloc(n*n2*LEN);
	tt=(double*)malloc(ndf*n2*LEN);
	ss1=(double*)malloc(n*np*LEN);
	ss2=(double*)malloc(n*np*np*LEN);
	eta=(double*)malloc(n*np*LEN);
	q=(double*)malloc(n*np*LEN);
	vt=(double*)malloc(np*np*LEN);
	xi=(double*)malloc(nc*np*LEN);


	//relates the matrices with the vector parameters
	if (n1 > 0)
        for (i = 0; i < n; i++)
            for (j = 0; j < n1; j++)
                *(a+i*n1+j)= x[i + n * j];

    if (n2 > 0)
        for (i = 0; i < n; i++)
            for (j = 0; j < n2; j++)
                *(aa+i*n2+j)= x2[i + n * j];
    if (n2 > 0)
        for (i = 0; i < ndf; i++)
            for (j = 0; j < n2; j++)
                *(tt+i*n2+j)= ft[i + ndf * j];

    // initialization
	for (i = 0; i < n; i++)
	{
		pi[i] = 0;
		for (j = 0; j < np; j++)
			*(eta+i*np+j) = *(q+i*np+j) = 0;
	}

	for (i = 0; i < n; i++)
    {
        ss0[i] = 0;
        for (j = 0; j < np; j++)
        {
            *(ss1+i*np+j)= 0;
            for (k = 0; k < np; k++)
                *(ss2+i*np*np+j*np+k) = 0;
        }
    }

	for (i = 0; i < np; i++)
        for (j = 0; j < np; j++)
			*(vt+i*np+j) = 0;

	for (i = 0; i < nc; i++)
        for (j = 0; j < np; j++)
			*(xi+i*np+j) = 0;


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
                    z[j1] = *(a+j*n1+j1);
                    zb += b[j1] * (*(a+j*n1+j1));
                }

            if (n2 > 0)
                for (k = 0; k < n2; k++)
                {
                    z[n1+k] = (*(aa+j*n2+k)) * (*(tt+count*n2+k));
                    zb += b[n1+k] * (*(aa+j*n2+k)) * (*(tt+count*n2+k));
                }

            if (t2[j] >= t2[i])
                wye = exp(zb);
            else
                wye = exp(zb) * wt[i] / wt[j];

            ss0[i] += wye;
            for (k = 0; k < np; k++)
            {
                *(ss1+i*np+k) += wye * z[k];
                for (j1 = 0; j1 < np; j1++)
                    *(ss2+i*np*np+k*np+j1) += wye * z[k] * z[j1];
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
            *(eta+i*np+k) = 0;

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
                    z[j1] = *(a+i*n1+j1);
                    zb += b[j1] * (*(a+i*n1+j1));
                }

            if (n2 > 0)
                for (k = 0; k < n2; k++)
                {
                    z[n1+k] = (*(aa+i*n2+k))*(*(tt+count1*n2+k));
                    zb += b[n1+k] * (*(aa+i*n2+k))*(*(tt+count1*n2+k));
                }


            if (t2[j] <= t2[i])
                wye = exp(zb);
            else if (t2[j] > t2[i] && ici[i]>1)
                wye = exp(zb) * wt[j] / wt[i];
            else continue;

            for (k = 0; k < np; k++)
                *(eta+i*np+k) -= (z[k] - (*(ss1+j*np+k))/ss0[j]) * wye / ss0[j];

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
					z[j] = *(a+i*n1+j);
					zb += b[j] * (*(a+i*n1+j));
				}
			if (n2 > 0)
				for (k = 0; k < n2; k++)
				{
					z[n1+k]= (*(aa+i*n2+k))*(*(tt+count*n2+k));
					zb += b[n1+k] * (*(aa+i*n2+k))*(*(tt+count*n2+k));
				}

			for (k = 0; k < np; k++)
				*(eta+i*np+k) += z[k] - (*(ss1+i*np+k))/ss0[i];

			//information
			for (j = 0; j < np; j++)
				for (k = 0; k < np; k++)
					*(vt+j*np+k) += (*(ss2+i*np*np+j*np+k))/ss0[i] - *(ss1+i*np+j) * (*(ss1+i*np+k))/ss0[i]/ss0[i];

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
							z[j] = *(a+j2*n1+j);
							zb += b[j] *(*(a+j2*n1+j));
						}
					if (n2 > 0)
						for (k = 0; k < n2; k++)
						{
							z[n1+k]= (*(aa+j2*n2+k)) * (*(tt+(count2-1)*n2+k));
							zb += b[n1+k] * (*(aa+j2*n2+k)) * (*(tt+(count2-1)*n2+k));
						}

					wye += exp(zb)* wt[j1] / wt[j2];
					for (k = 0; k < np; k++)
						wyez[k] += exp(zb)* wt[j1] / wt[j2] * z[k];
				}

				for (k = 0; k < np; k++)
					*(q+i*np+k) += (wyez[k] - (*(ss1+j1*np+k))/ss0[j1] * wye) /ss0[j1];
			}
		}  
	}

   //////////////////////////////////////////
   //summarize the variance and output results
	for (i = 0; i < n; i++)
        for (k = 0; k < np; k++)
		{
			if (ici[i]==0) *(eta+i*np+k) += (*(q+i*np+k))/pi[i];
		    for (j = 0; j < n; j++)
			{
				    if (t2[j] > t2[i]) break;
				    if (ici[j] == 0)  *(eta+i*np+k) -= (*(q+j*np+k))/pi[j]/pi[j];
			 }
		}


    for (j = 0; j < np; j++)
        for (k= 0; k < np; k++)
            v[j  + k * np] = *(vt+j*np+k) ;   

	for (i=0; i< np*np; i++)
		v2[i]=0;

	if (nc == n)
	{
		for (j = 0; j < np; j++)
			for (k= 0; k < np; k++)
				for (i = 0; i < n; i++)
					v2[j  + k * np] += (*(eta+i*np+j)) * (*(eta+i*np+k));
	}
	else
	{
		for (i = 0; i < nc; i++)
			for (j = 0; j < n; j++)
			{
				if (cluster[j] != i+1) continue;
				for (k = 0; k < np; k++)
					*(xi+i*np+k) += *(eta+j*np+k);
			}

		for (j = 0; j < np; j++)
			for (k= 0; k < np; k++)
				for (i = 0; i < nc; i++)
					v2[j  + k * np] += (*(xi+i*np+j)) * (*(xi+i*np+k));
	}


	free(a);free(aa);free(tt);free(vt);free(ss1);free(ss2);free(eta);free(q);free(xi);
}
