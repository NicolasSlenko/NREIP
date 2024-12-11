function [cl, cd, converge]=xfoil(X,Y,ADEG,RE,MACH)

% ////--------------------------------------
% //     xfoil function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

     
     p=initialize();
     % Newton iteration
     p.itmax=30;
	%//---- set minf, reinf, based on current cl-dependence
	[p,minf_cl, reinf_cl]=mrcl(1.0, p);
	
	%//---- set various compressibility parameters from minf
	p=comset(p);
	
% 	return true;
% }
p.xb=X;
p.yb=Y;
p.nb=length(X);



      p.area = 0.0;
      for I=1: p.nb
        IP = I+1;
        if (I==p.nb) 
            IP = 1;
        end
        p.area = p.area + 0.5*(p.yb(I)+p.yb(IP))*(p.xb(I)-p.xb(IP));
      end
% C
      if(p.area>=0.0) 
       p.lclock = false;
%        WRITE(*,1010) NB
      else
% C----- if area is negative (clockwise order), reverse coordinate order
       p.lclock = true;
%        WRITE(*,1011) NB
       for I=1: p.nb/2
         p.xtmp = p.xb(p.nb-I+1);
         p.ytmp = p.yb(p.nb-I+1);
         p.xb(p.nb-I+1) = p.xb(I);
         p.yb(p.nb-I+1) = p.yb(I);
         p.xb(I) = p.xtmp;
         p.yb(I) = p.ytmp;
%    55  CONTINUE
%       ENDIF
       end
      end
      
      p.sb= scalc(p.xb,p.yb,p.nb);
      p.xbp= segspl(p.xb,p.sb,p.nb);
      p.ybp= segspl(p.yb,p.sb,p.nb);
      
	  [p.w1,p.sble,p.chordb,p.areab,...
		p.radble,p.angbte,p.ei11ba,p.ei22ba,p.apx1ba,p.apx2ba,...
		p.ei11bt,p.ei22bt,p.apx1bt,p.apx2bt]=geopar(p.xb,p.xbp,p.yb,p.ybp,p.sb,p.nb);

      
      
      p.xble = seval(p.sble,p.xb,p.xbp,p.sb,p.nb);
      p.yble = seval(p.sble,p.yb,p.ybp,p.sb,p.nb);
      p.xbte = 0.5*(p.xb(1) + p.xb(p.nb));
      p.ybte = 0.5*(p.yb(1) + p.yb(p.nb));

% C---- wipe out old flap hinge location
      p.xbf = 0.0;
      p.ybf = 0.0;
      p.lbflap = false;
      
      [p]=abcopy(p);
   
   
         [p.imax,p.amax]= cang(p.x,p.y,p.n);
      
      
         p.lalfa = true;
         p.alfa = ADEG*pi/180;
         p.qinf = 1.0;
         
         
         [p]=specal(p);
         	CL=p.cl;
    CD=p.cdp;
    CM=p.cm;
                    if (abs(p.alfa-p.awake) > 1.0e-5) 
           p.lwake  = false;
         end
         if (abs(p.alfa-p.avisc) > 1.0e-5)
         p.lvconv = false;
         end
%          if (abs(p.minf-p.mvisc) > 1.0e-5) 
%              p.lvconv = false;
%          end
        if (RE==0 && MACH==0)
         p.lvisc=false;
        else 
             p.lvisc=true;
        end
         p.reinf1=RE;
         p.minf1=MACH;
         p.retyp=1;
         p.matyp=1;
         
      if (p.lvisc==true) 
%          [CL]= VISCAL(ITMAX,LWAKE,N,S,WAKLEN,CHORD,X,Y,XP,YP,NX,NY,SHARP,ANTE,ASTE,DSTE,APANEL,SIG,GAMU,GAM,ALFA,QINVU,LIPAN,LBLINI,LVCONV,BIJ,AIJ,AIJPIV,QINF,MINF,LADIJ,AIJ,AIJPIV,BIJ,LALFA,CL,RETYP,MATYP,MINF1,REINF1,GAMM1,GAMMA,ACRIT,IDAMP,XSTRIP,XLE,YLE,SLE);
        [p, converge]=viscal(p.itmax,p);
        disp(converge);
%          CL=p.cl;
%          CD=p.cd;
%          CM=p.cm;
      end
%               p.nrp = 21;
%         WRITE(*,*) 'Using default number of profiles:', NPR
%        ENDIF
% C
%        IF(NPR.GT.1) THEN
% % C------ set NPR points along surface, offset slightly for the locating logic
%         p.doff = 0.00001*(p.s(p.n)-p.s(1));
%         for ipr = 1: p.nrp
%           frac = (ipr-1)/(p.nrp-1);
%           spr = p.s(1) + (p.s(p.n)-p.s(1))*frac;
%           p.xpr(ipr) = seval(spr,p.x,p.xp,p.s,p.n) + p.doff*deval(spr,p.y,p.yp,p.s,p.n);
%           p.ypr(ipr) = seval(spr,p.y,p.yp,p.s,p.n) - p.doff*deval(spr,p.x,p.xp,p.s,p.n);
%         end
%         [p]=dplot(p);
       
    cl = p.cl;
    cd = p.cd;

end
function [p]=dplot(p)
% ////--------------------------------------
% //     dplot function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

%       INCLUDE 'XFOIL.INC'     
% C-----------------------------------------------------------
% C     Plots analytical profiles at specified points.
% C     If p.nrp=0, then cursor-selected points are requested.
% C-----------------------------------------------------------
%       DIMENSION p.xpr(*), p.ypr(*)
% C
%       CHARACTER*1 KCHAR
%       LOGICAL lcrs, p.turb
%       LOGICAL LGUI
% C
%       CALL GETCOLOR(ICOL0)
% C
%       lcrs = p.nrp .LE. 0
% C
%       if (lcrs) 
%        kdone  = 1
%        xdwin = p.xpage - 2.0*p.xmarg
%        ydwin = p.ypage - 2.0*YMARG
%        X1 = p.xmarg + 0.91*xdwin
%        X2 = p.xmarg + 0.99*xdwin
%        Y1 = YMARG + 0.01*ydwin
%        Y2 = YMARG + 0.05*ydwin
%        CALL NEWPEN(5)
%        CALL GUIBOX(kdone, X1,X2,Y1,Y2, 'GREEN'   , ' Done ')
% C
%        WRITE(*,*) ' '
%        WRITE(*,*) 'Locate profiles with cursor, type "D" when done...'
%        p.nrp = 12345
% C
%       ELSE
%        p.nrp = p.nrp
% C
%       ENDIF
% C
% C---- go over profiles ...
p.uprwt=0.02;
      for IPR=1: p.nrp
% C
%         IF(lcrs) THEN 
% C------- get cursor plot coordinates
%          CALL GETCURSORXY(p.xc,p.yc,KCHAR)
%          IF(INDEX('Dd',KCHAR).NE.0 .OR. LGUI(kdone,p.xc,p.yc)) THEN
%           RETURN
%          ENDIF
% C
% C------- transform to airfoil coordinates
%          p.xc = p.xc/p.faca - p.xofa
%          p.yc = p.yc/p.faca - p.yofa
% C
%         ELSE
         p.xc = p.xpr(IPR);
         p.yc = p.ypr(IPR);
% C
%         ENDIF
% C
% C------ find nearest airfoil surface point
        p.rsqmin = 1.0E23;
        p.ismin = 0;
        p.iblmin = 0;
%         p.doff = 0.00001*(p.s(N)-p.s(1))
        for is= 1: 2
          for ibl= 2: p.iblte(is)
            I = p.ipan(ibl,is);
            p.xsurf = p.x(I) + p.doff*p.yp(I);
            p.ysurf = p.y(I) - p.doff*p.xp(I);
            p.rsq = (p.xc-p.xsurf)^2 + (p.yc-p.ysurf)^2;
            if (p.rsq <= p.rsqmin) 
             p.rsqmin = p.rsq;
             p.ismin = is;
             p.iblmin = ibl;
            end
          end
        end
% C
        is = p.ismin;
        ibl = p.iblmin;
% C
        I = p.ipan(ibl,is);
        p.crsp = (p.xc-p.x(I))*p.ny(I) - (p.yc-p.y(I))*p.nx(I);
        if(is==2) 
            p.crsp = -p.crsp;
        end
        if (p.crsp>0.0) 
         p.iblp = ibl+1;
         p.iblo = ibl;
        else
         p.iblp = ibl;
         p.iblo = ibl-1;
        end
        ISP = is;
        ISO = is;
% C
        if (p.iblp>p.iblte(is)) 
         p.iblp = p.iblte(is);
         p.iblo = p.iblp-1;
         ibl = p.iblte(is);
        elseif(p.iblo<2) 
         p.iblo = 2;
         if (ISO==1) 
          ISO = 2;
         else
          ISO = 1;
         end
        end
% C
        IP = p.ipan(p.iblp,ISP);
        IO = p.ipan(p.iblo,ISO);
% C
% C------ set interpolation fraction at profile location
        dx = p.x(IP) - p.x(IO);
        dy = p.y(IP) - p.y(IO);
        vx = p.xc - p.x(IO);
        vy = p.yc - p.y(IO);
        frac = (dx*vx + dy*vy)/(dx*dx+dy*dy);
        frac = min( max( frac , 0.0 ) , 1.0 );
% C
% C------ set averaged displacement vector at profile location
        p.ca = frac*p.ny(IP) + (1.0-frac)*p.ny(IO);
        p.sa = frac*p.nx(IP) + (1.0-frac)*p.nx(IO);
        p.csmod = sqrt(p.ca^2 + p.sa^2);
        p.ca = p.ca/p.csmod;
        p.sa = p.sa/p.csmod;
% C
        p.xo = frac*p.x(IP) + (1.0-frac)*p.x(IO);
        p.yo = frac*p.y(IP) + (1.0-frac)*p.y(IO);
% C
        DS = frac*p.dstr(p.iblp,ISP) + (1.0-frac)*p.dstr(p.iblo,ISO);
        p.th = frac*p.thet(p.iblp,ISP) + (1.0-frac)*p.thet(p.iblo,ISO);
        UE = frac*p.uedg(p.iblp,ISP) + (1.0-frac)*p.uedg(p.iblo,ISO);
% C
        XI = frac*p.xssi(p.iblp,ISP) + (1.0-frac)*p.xssi(p.iblo,ISO);
        p.turb = XI > p.xssitr(is);
% C
% C------ 1 / (total enthalpy)
        p.hstinv = p.gamm1*(p.minf/p.qinf)^2 / (1.0 + 0.5*p.gamm1*p.minf^2);
% C
% C------ Sutherland's const./To   (assumes stagnation conditions are at STP)
%         p.hvrat = 0.35
% C
% C------ fill Rtheta arrays
        UEC = UE * (1.0-p.tklam) / (1.0 - p.tklam*(UE/p.qinf)^2);
        HERAT = (1.0 - 0.5*p.hstinv*UEC ^2)...
            / (1.0 - 0.5*p.hstinv*p.qinf^2);
        RHOE = HERAT ^ (1.0/p.gamm1);
        AMUE = sqrt(HERAT^3) * (1.0+p.hvrat)/(HERAT+p.hvrat);
        p.rtheta = p.reinf * RHOE*UE*p.th/AMUE;
% C
        p.amsq = UEC*UEC*p.hstinv / (p.gamm1*(1.0 - 0.5*UEC*UEC*p.hstinv));
        [p.hk, DUMMY, DUMMY]=hkin( DS/p.th, p.amsq);
% C
%         WRITE(*,9100) p.xo,p.yo, DS, p.rtheta, p.hk
%  9100   FORMAT(1X,'x y =', 2F8.4,'    Delta* =', G12.4,
%      &         '    Rtheta =', F10.2,'    Hk =', F9.4)
% C
        if (is==1) 
         p.udir = 1.0;
        else
         p.udir = -1.0;
        end
% C
        p.uei = UE/p.qinf;
        p.un = 0.0;
%         CALL NEWCOLORNAME('green')
        
        p.uprwts = p.uprwt*0.5*(p.s(p.n)-p.s(1));
%         CALL PRPLOT(p.xo,p.yo,p.th,p.uei,p.un,p.hk,p.rtheta,p.amsq,p.turb,
%      &              -XOFA,-YOFA,FACA,p.uprwts,p.sa,p.ca,p.udir)
%    50 CONTINUE
      end
% C
%       CALL NEWCOLOR(ICOL0)
%       CALL PLFLUSH
% C
%       RETURN
%      END ! DPLOT
end
function [p, converge]=viscal(NITER,p)
% ////--------------------------------------
% //     viscal function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% ////--------------------------------------
% //     converges viscous operating point
% ////--------------------------------------
% 	int ibl;
% 
% //---- calculate wake trajectory from current inviscid solution if necessary
	if (p.lwake==false)
        [p]=xyWake(p);
    end
	
% //	---- set velocities on wake from airfoil vorticity for alpha=0, 90
	[p]=qwcalc(p);

% //	---- set velocities on airfoil and wake for initial alpha
	[p]=qiset(p);

	if(p.lipan==false)% {
		
		if(p.lblini) 
            p=gamqv(p);
        end
		
% 		//	----- locate stagnation point arc length position and panel index
		[p]=stfind(p);
		
% 		//	----- set  bl position -> panel position  pointers
		[p]=iblpan(p);
		
% 		//	----- calculate surface arc length array for current stagnation point location
		[p]=xicalc(p);
		
% 		//	----- set  bl position -> system line  pointers
		[p]=iblsys(p);

    end
	
% //	---- set inviscid bl edge velocity uinv from qinv
	[p]= uicalc(p);
	
	if p.lblini==false% {
% 		//	----- set initial ue from inviscid ue
		for ibl=1: p.nbl(1)%{
			p.uedg(ibl,1) = p.uinv(ibl,1);
        end
		for ibl=1:p.nbl(2)
			p.uedg(ibl,2) = p.uinv(ibl,2);
        end
    end
	
	if (p.lvconv )% {
% 		//	----- set correct cl if converged point exists
		p=qvfue(p);
% /*
% 		if(lvisc)%{
% 			if(!cpcalc(n+nw,qvis,qinf,minf,cpv)){
% 				return false;
% 			}
% 			if(!cpcalc(n+nw,qinv,qinf,minf,cpi)){
% 				return false;
% 			}
% 		}
% 		else if(!cpcalc(n,qinv,qinf,minf,cpi)){
% 			return false;
% 		}
% */
		if(p.lvisc)% {
			p.cpv=cpcalc(p.n+p.nw,p.qvis,p.qinf,p.minf);
			p.cpi=cpcalc(p.n+p.nw,p.qinv,p.qinf,p.minf);
		
        else
            p.cpi=cpcalc(p.n,p.qinv,p.qinf,p.minf);
        end
		p=gamqv(p);
		p=clcalc(p.xcmref,p.ycmref,p);
		p=cdcalc(p);
    end
	
% //	---- set up source influence matrix if it doesn't exist
	if(p.lwdij==false || p.ladij==false) 
        [p]=qdcalc(p);
    end

    eps1  = 1; %=0.005;
    for ITER=1: NITER
        
        [p]=ViscousIter(p);
            if(p.rmsbl < eps1) 
		        p.lvconv = true;
		        p.avisc = p.alfa;	
		        p.mvisc = p.minf;
                p.cpi=cpcalc(p.n+p.nw,p.qinv,p.qinf,p.minf);
	            p.cpv=cpcalc(p.n+p.nw,p.qvis,p.qinf,p.minf);
	                if(p.lflap)
                        p=mhinge(p)
                    end
		        disp('----------CONVERGED----------\r\n\r\n')
                converge = 1;
                return
            end
    end
   disp( 'VISCAL:  Convergence failed')
   converge = 0;

end
function [p]=ViscousIter(p)
% ////--------------------------------------
% //     ViscousIter function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------


% //	Performs one iteration
	
% 	CString str;

    
	[p]=setbl(p);
    %//	------ fill newton system for bl variables
    
	[p]=blsolve(p);%//	------ solve newton system with custom solver	

	p=update(p);%//	------ update bl variables


	if(p.lalfa) %{//	------- set new freestream mach, re from new cl
        [p,p.minf_cl,p.reinf_cl]=mrcl(p.cl,p);
% 		mrcl(cl, minf_cl, reinf_cl);
% 		comset();
        [p]=comset(p);
% 	}%
    else%{//	------- set new inviscid speeds qinv and uinv for new alpha
		p=qiset(p);
		p=uicalc(p);
    end
% 	}

	
	p=qvfue(p);%//	------ calculate edge velocities qvis(.) from uedg(..)
	p=gamqv(p);%//	------ set gam distribution from qvis
	p=stmove(p);%//	------ relocate stagnation point
	
% 	//	------ set updated cl,cd
[p]=clcalc(p.xcmref, p.ycmref,p);
% 	clcalc(xcmref,ycmref);
	p=cdcalc(p);
	
% 	//	------ display changes and test for convergence
% 	if(rlx<1.0) str.Format("     rms:%.2e   max:%.2e at %d %d   rlx:%.3f\r\n",rmsbl, rmxbl, imxbl,ismxbl,rlx);
% 	else if(abs(rlx-1.0)<0.001) str.Format("     rms:%.2e  max:%.2e at %d %d \r\n", rmsbl, rmxbl,imxbl,ismxbl);
% 	if(m_bTrace)pXFile->WriteString(str);

	p.cdp = p.cd - p.cdf;

% 	str.Format("     a=%.3f    cl=%.4f \r\n     cm=%.4f  cd=%.5f => cdf=%.5f cdp=%.5f\r\n\r\n", alfa/dtor, cl, cm, cd, cdf, cdp);
% 	if(m_bTrace)pXFile->WriteString(str);

% 	int pos = str.Find("QN",0);
end
function cp=cpcalc( n,  q,  qinf,  minf)
% ////--------------------------------------
% //     cpcalc function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% {
% 	//---------------------------------------------
% 	//     sets compressible cp from speed.
% 	//---------------------------------------------
% 
% 	bool denneg;
% 	 cpinc, den;
	 beta = sqrt(1.0 - minf*minf);
	 bfac = 0.5*minf*minf / (1.0 + beta);
	
	denneg = false;
	
	for (i=1:n)
% 	{
		cpinc = 1.0 - (q(i)/qinf)*(q(i)/qinf);
		den = beta + bfac*cpinc;
		cp(i) = cpinc / den;
		if(den <= 0.0) 
            denneg = true;
        end
    end
	
	if(denneg)  
% 	{
% 		CString str;
		disp(' CpCalc: local speed too large \r\n Compressibility corrections invalid \r\n');
% //		AfxMessageBox(str);
% 		if(pXFile) pXFile->WriteString(str);
		return;
    end
	
% 	return true;
end
function p=xyWake(p)
% ////--------------------------------------
% //     xyWake function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% {
% 
% //-----------------------------------------------------
% //     sets wake coordinate array for current surface 
% //     vorticity and/or mass source distributions.
% //-----------------------------------------------------
% 	double ds, ds1, p.sx, p.sy, p.smod;
% 	double psi, psi_x,psi_y;
% //
% 	//TRACE("calculating wake trajectory ...\p.n");
% 		CString str;
		disp(' Calculating wake trajectory ...\r\p.n');
% 		if(m_bTrace)pXFile->WriteString(str);
% 	//
% 	//--- number of wake points
	p.nw = p.n/8 + 2;
% 	if(p.nw>IWX) {
% //		AfxMessageBox(" xywake: array size (iwx) too small.  last wake point index reduced.", MB_ICONSTOP | MB_OK);
% 		CString str;
% 		str.Format(" XYWake: array size (IWX) too small.\p.n  Last wake point index reduced.");
% 		if(m_bTrace)pXFile->WriteString(str);
% 		p.nw = IWX;
% 	}
	
	ds1 = 0.5*(p.s(2) - p.s(1) + p.s(p.n) - p.s(p.n-1));
	p.snew=setexp(ds1,p.waklen*p.chord,p.nw);
	p.snew=[ones(1,p.n) p.snew];
	p.xte = 0.5*(p.x(1)+p.x(p.n));
	p.yte = 0.5*(p.y(1)+p.y(p.n));
	
% 	//-- set first wake point a tiny distance behind te
	 i = p.n+1;
	p.sx = 0.5*(p.yp(p.n) - p.yp(1));
	p.sy = 0.5*(p.xp(1) - p.xp(p.n));
	p.smod = sqrt(p.sx*p.sx + p.sy*p.sy);
	p.nx(i) = p.sx / p.smod;
	p.ny(i) = p.sy / p.smod;
	p.x(i) = p.xte - 0.0001*p.ny(i);
	p.y(i) = p.yte + 0.0001*p.nx(i);
	p.s(i) = p.s(p.n);
	
% 	//---- calculate streamfunction gradient components at first point
	[psi,psi_x,p]=psilin(i,p.x(i),p.y(i),1.0,0.0,false,false,p);
	[psi,psi_y,p]=psilin(i,p.x(i),p.y(i),0.0,1.0,false,false,p);

% 	//---- set unit vector normal to wake at first point
	p.nx(i+1) = -psi_x / sqrt(psi_x*psi_x + psi_y*psi_y);
	p.ny(i+1) = -psi_y / sqrt(psi_x*psi_x + psi_y*psi_y);
	
% 	//---- set angle of wake panel normal
	p.apanel(i) = atan2( psi_y , psi_x );
	
% 	//---- set rest of wake points
	for(i=p.n+2:p.n+p.nw)
		ds = p.snew(i) - p.snew(i-1);
		
% 		//------ set new point ds downstream of last point
		p.x(i) = p.x(i-1) - ds*p.ny(i);
		p.y(i) = p.y(i-1) + ds*p.nx(i);
		p.s(i) = p.s(i-1) + ds;
		
		if(i~=p.n+p.nw) 
			
% 			//------- calculate normal vector for next point
			[psi,psi_x,p]=psilin(i,p.x(i),p.y(i),1.0,0.0,false,false,p);
			[psi,psi_y,p]=psilin(i,p.x(i),p.y(i),0.0,1.0,false,false,p);
			
			p.nx(i+1) = -psi_x / sqrt(psi_x*psi_x + psi_y*psi_y);
			p.ny(i+1) = -psi_y / sqrt(psi_x*psi_x + psi_y*psi_y);
			
% 			//------- set angle of wake panel normal
			p.apanel(i) = atan2( psi_y , psi_x );
        end
    end
	
% 	//---- set wake presence flag and corresponding alpha
	p.lwake = true;
	p.awake =  p.alfa;
	
% 	//---- old source influence matrix is invalid for the new wake geometry
	p.lwdij = false;
	
end
function s=setexp( ds1,  smax, nn)
% ////--------------------------------------
% //     setexp function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% 	//........................................................
% 	//     sets geometriy stretched array s:
% 	//
% 	//       s(i+1) - s(i)  =  r * [s(i) - s(i-1)]
% 	//
% 	//       s     (output)  array to be set  
% 	//       ds1   (input)   first s increment:  s[2] - s[1]
% 	//       smax  (input)   final s value:      s(nn)
% 	//       nn    (input)   number of points
% 	//........................................................
% 	CString str;
% 
% 	int nex, iter;
% 	 sigma, rnex, rni, aaa, bbb, ccc;
% 	 disc, ratio, sigman, res;
% 	 dresdr, dratio, ds;
% 	
goto=0;
	sigma = smax/ds1;
	nex = nn-1;
	rnex = nex;
	rni = 1.0/rnex;
	
% 	//-- solve quadratic for initial geometric ratio guess
	aaa = rnex*(rnex-1.0)*(rnex-2.0) / 6.0;
	bbb = rnex*(rnex-1.0) / 2.0;
	ccc = rnex - sigma;
	
	disc = bbb*bbb - 4.0*aaa*ccc;
	disc = max(0.0, disc);
	
% 	if(nex<=1) {
% //		AfxMessageBox("setexp: cannot fill array.  n too small", MB_ICONSTOP | MB_OK);
% 		CString str;
% 		str.Format("setexp: cannot fill array.  n too small\r\n");
% 		if(m_bTrace)pXFile->WriteString(str);
% 		return false;
% 	}
% 	else {
		if(nex==2)    
            ratio = -ccc/bbb  +  1.0;
        else
            ratio = (-bbb + sqrt(disc))/(2.0*aaa)  +  1.0;
        end
	if(ratio==1.0) 
        goto=11;
    end
		
% 	//-- newton iteration for actual geometric ratio
  if goto~=11
	for (iter=1:100)
		sigman = ((ratio^nex) - 1.0) / (ratio - 1.0);
		res = (sigman^rni) - (sigma^rni);
		dresdr = rni*(sigman^rni)...
			* (rnex*(ratio^(nex-1)) - sigman) ...
			/ ((ratio^nex) - 1.0);
		
		dratio = -res/dresdr;
		ratio = ratio + dratio;
		
		if(abs(dratio) < 1.0e-5) 	
            goto=11;
            break
        end
    end
  end
  if goto~=11
	disp('setexp: convergence failed.  continuing anyway ...');
% 	str.Format("Setexp: Convergence failed.  Continuing anyway ...\r\n");
% 	if(m_bTrace)pXFile->WriteString(str);
  end


% 	//-- set up stretched array using converged geometric ratio
% stop11:
	s(1) = 0.0;
	ds = ds1;
	for (n=2: nn)
		s(n) = s(n-1) + ds;
		ds = ds*ratio;
    end
% 	return true;
end
function p=qwcalc(p)
% ////--------------------------------------
% //     qwcalc function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% {
% 
% 	//---------------------------------------------------------------
% 	//     sets inviscid tangential velocity for alpha = 0, 90
% 	//     on wake due to freestream and airfoil surface vorticity.
% 	//---------------------------------------------------------------
% 	double psi, psi_ni;
% 	
% 	//---- first wake point (same as te)
	p.qinvu(p.n+1,1) = p.qinvu(p.n,1);
	p.qinvu(p.n+1,2) = p.qinvu(p.n,2);
	
% 	//---- rest of wake
	for (i=p.n+2:p.n+p.nw)
		[psi,psi_ni,p]=psilin(i,p.x(i),p.y(i),p.nx(i),p.ny(i),false,false,p);
		p.qinvu(i,1) = p.qtan1;
		p.qinvu(i,2) = p.qtan2;
    end
	
% 	return true;
end
function p=stfind(p)
% ////--------------------------------------
% //     stfind function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% {
% 	//-----------------------------------------
% 	//     locates stagnation point arc length 
% 	//     location p.sst and panel index ist.
% 	//-----------------------------------------
% 	double p.dgam, p.ds;
% 	int i;

	p.bFound = false;
	for(i=1:p.n-1)
		if(p.gam(i)>=0.0 && p.gam(i+1)<0.0) 
			p.bFound = true;
			break;
        end
    end

	if(p.bFound==false)
		disp('stfind: stagnation point not found. continuing ...');
% 		CString str;
% 		str.Format("stfind: Stagnation point not found. Continuing ...\r\n");
% 		if(m_bTrace)pXFile->WriteString(str);
		i = p.n/2;
    end

% //stop11:	
	p.ist = i;
	p.dgam = p.gam(i+1) - p.gam(i);
	p.ds = p.s(i+1) - p.s(i);
	
% 	//---- evaluate so as to minimize roundoff for very small p.gam(i) or p.gam(i+1)
	if(p.gam(i) < -p.gam(i+1))
		p.sst = p.s(i)   - p.ds*(p.gam(i)  /p.dgam);
	else
		p.sst = p.s(i+1) - p.ds*(p.gam(i+1)/p.dgam);
    end

% 	//---- tweak stagnation point if it falls right on a node (very unlikely)
	if(p.sst <= p.s(i)  ) 
		p.sst = p.s(i)   + 0.0000001;
    end
	if(p.sst >= p.s(i+1)) 
		p.sst = p.s(i+1) - 0.0000001;
    end
	
	p.sst_go = (p.sst  - p.s(i+1))/p.dgam;
	p.sst_gp = (p.s(i) - p.sst   )/p.dgam;
	
% 	return true;
% 
end
function p=iblpan(p)
% ////--------------------------------------
% //     iblpan function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% {
% 	//-----------------------------------------------------------
% 	//     sets  bl location -> panel location  pointer array p.ipan
% 	//-----------------------------------------------------------
% 	int iblmax,is, ibl, i, iw;
% 
% 	//-- top surface first
	is = 1;
	
	ibl = 1;
	for(i=p.ist:-1: 1)
		ibl = ibl+1;
		p.ipan(ibl,is) = i;
		p.vti(ibl,is) = 1.0;
    end
	
	p.iblte(is) = ibl;
	p.nbl(is) = ibl;
	
% 	//-- bottom surface next
	is = 2;
	ibl = 1;
	for(i=p.ist+1:p.n)
		ibl = ibl+1;
		p.ipan(ibl,is) = i;
		p.vti(ibl,is) = -1.0;
    end
	
% 	//-- wake
	p.iblte(is) = ibl;
	
	for(iw=1:p.nw)
		
		i = p.n+iw;
		ibl = p.iblte(is)+iw;
		p.ipan(ibl,is) = i;
		p.vti(ibl,is) = -1.0;
    end
	
	p.nbl(is) = p.iblte(is) + p.nw;
	
% 	//-- upper wake pointers (for plotting only)
	for(iw=1:p.nw)
		p.ipan(p.iblte(1)+iw,1) = p.ipan(p.iblte(2)+iw,2);
		p.vti(p.iblte(1)+iw,1) = 1.0;
    end
	p.iblmax = max(p.iblte(1),p.iblte(2)) + p.nw;
% 	if(iblmax>IVX) {
% //		AfxMessageBox("iblpan :  ***  bl array overflow", MB_ICONSTOP | MB_OK);
% 		CString str;
% 		str.Format("iblpan :  ***  bl array overflow");
% 		if(m_bTrace)pXFile->WriteString(str);
% 
% 		str.Format("Increase IVX to at least %d\r\n", iblmax);
% 		if(m_bTrace)pXFile->WriteString(str);
% 		return false;
% 	}
	
	p.lipan = true;
end
function p=xicalc(p)
% ////--------------------------------------
% //     xicalc function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% % 	//-------------------------------------------------------------
% % 	//     sets bl arc length array on each airfoil side and wake
% % 	//-------------------------------------------------------------
% % 	double p.telrat, p.crosp, p.dwdxte, aa, bb, zn;
% 	int i, ibl;
	is = 1;
	
	p.xssi(1,is) = 0.0;
	
	for (ibl=2:p.iblte(is))
		i = p.ipan(ibl,is);
		p.xssi(ibl,is) = p.sst - p.s(i);
    end
	
	is = 2;
	
	p.xssi(1,is) = 0.0;
	
	for (ibl=2:p.iblte(is))
		
		i = p.ipan(ibl,is);
		p.xssi(ibl,is) = p.s(i) - p.sst;
    end
	
	ibl = p.iblte(is) + 1;
	p.xssi(ibl,is) = p.xssi(ibl-1,is);
	
	for (ibl=p.iblte(is)+2:p.nbl(is))
		
		i = p.ipan(ibl,is);
		p.xssi(ibl,is) = p.xssi(ibl-1,is)...
			+ sqrt((p.x(i)-p.x(i-1))* (p.x(i)-p.x(i-1))...
			+ (p.y(i)-p.y(i-1))*(p.y(i)-p.y(i-1)));
    end
	
% 	//---- trailing edge flap length to te gap ratio
	p.telrat = 2.50;
	
% 	//---- set up parameters for te flap cubics
	
% 	//   p.dwdxte = p.yp(1)/p.xp(1) + p.yp(p.n)/p.xp(p.n)    !!! bug  2/2/95
	
	p.crosp = (p.xp(1)*p.yp(p.n) - p.yp(1)*p.xp(p.n))...
		/ sqrt(  (p.xp(1)*p.xp(1) + p.yp(1)*p.yp(1))...
		*(p.xp(p.n)*p.xp(p.n) + p.yp(p.n)*p.yp(p.n)) );
	p.dwdxte = p.crosp / sqrt(1.0 - p.crosp*p.crosp);
	
% 	//---- limit cubic to avoid absurd te gap widths
	p.dwdxte = max(p.dwdxte,-3.0/p.telrat);
	p.dwdxte = min(p.dwdxte, 3.0/p.telrat);
	
	aa =  3.0 + p.telrat*p.dwdxte;
	bb = -2.0 - p.telrat*p.dwdxte;
	
	if(p.sharp) 
		for (iw=1:p.nw)
			p.wgap(iw) = 0.0;
        end
  
	
    else
% 		//----- set te flap (wake gap) array
		is = 2;
		for (iw=1:p.nw)
			ibl = p.iblte(is) + iw;
			zn = 1.0 - (p.xssi(ibl,is)-p.xssi(p.iblte(is),is)) / (p.telrat*p.ante);
			p.wgap(iw) = 0.0;
			if(zn>=0.0)
                p.wgap(iw) = p.ante * (aa + bb*zn)*zn*zn;
            end
        end
    end
% 	return true;
end
function p=iblsys(p)
% ////--------------------------------------
% //     iblsys function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% {
% //---------------------------------------------
% //     sets the bl newton system line number
% //     corresponding to each bl station.
% //---------------------------------------------

	iv = 0;
	for (is=1:2)
		for (ibl=2:p.nbl(is))
			iv = iv+1;
			p.isys(ibl,is) = iv;
        end
    end
	
	p.nsys = iv;
% 	if(nsys>2*IVX) {
% //		AfxMessageBox("*** iblsys: bl system array overflow. ***", MB_ICONSTOP | MB_OK);
% 		CString str;
% 		str.Format("*** iblsys: bl system array overflow. ***");
% 		if(m_bTrace)pXFile->WriteString(str);
% 		return false;
% 	}
	
% 	return true;
end
function p=uicalc(p)
% ////--------------------------------------
% //     uicalc function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% 	//--------------------------------------------------------------
% 	//     sets inviscid ue from panel inviscid tangential velocity
% 	//--------------------------------------------------------------
	for (is=1:2)
        p.uinv  (1,is) = 0.0;
        p.uinv_a(1,is) = 0.0;
        for (ibl=2:p.nbl(is))
			i = p.ipan(ibl,is);
			p.uinv(ibl,is) = p.vti(ibl,is)*p.qinv  (i);
			p.uinv_a(ibl,is) = p.vti(ibl,is)*p.qinv_a(i);
        end
    end
	
% 	return true;
	
end
function p=qdcalc(p)
% ////--------------------------------------
% //     qdcalc function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% 	//-----------------------------------------------------
% 	//	   calculates source panel influence coefficient
% 	//	   matrix for current airfoil and wake geometry.
% 	//-----------------------------------------------------
% 	int i,j,k, iu;
% 	double psi, psi_n, sum;
% 	double bbb(IQX);
	
	disp('calculating source influence matrix ...');
% 	CString str;
% 	str.Format(" Calculating source influence matrix ...\r\n");
% 	if(m_bTrace)pXFile->WriteString(str);


	if(p.ladij==false) 
% 		//----- calculate source influence matrix for airfoil surface if it doesn't exist
		for (j=1:p.n)
			
% 			//------- multiply each dpsi/sig vector by inverse of factored dpsi/dgam matrix
% 			for (iu=0; iu<IQX; iu++) bbb(iu) = bij(iu,j);//arcds : create a dummy array
			p.bij(:,j)=baksub(p.n+1,p.aij,p.aijpiv,p.bij(:,j));
% 			for (iu=0; iu<IQX; iu++) bij(iu,j) = bbb(iu);
			
% 			//------- store resulting dgam/dsig = dqtan/dsig vector
			for (i=1:p.n)
				p.dij(i,j) = p.bij(i,j); 	
            end
        end
		p.ladij = true;
    end

% 	//---- set up coefficient matrix of dpsi/dm on airfoil surface
	for (i=1:p.n)
		[psi,psi_n,p]=pswlin(i,p.x(i),p.y(i),p.nx(i),p.ny(i),p);
		for (j=p.n+1:p.n+p.nw)
			p.bij(i,j) = -p.dzdm(j);
        end
    end
	
% 	//---- set up kutta condition (no direct source influence)
	for(j=p.n+1:p.n+p.nw) 
        p.bij(p.n+1,j) = 0.0;
    end
	
% 	//---- sharp te gamma extrapolation also has no source influence
	if(p.sharp) 
        for(j=p.n+1:p.n+p.nw)
            p.bij(p.n,j) = 0.0;
        end
    end
	
% 	//---- multiply by inverse of factored dpsi/dgam matrix
	for(j=p.n+1:p.n+p.nw)
% //		baksub(iqx,n+1,aijpiv,j);
% 		for (iu=0; iu<IQX; iu++) bbb(iu) = bij(iu,j);//arcds : create a dummy array
		
		p.bij(:,j)=baksub(p.n+1,p.aij,p.aijpiv,p.bij(:,j));
% 		for (iu=0; iu<IQX; iu++) bij(iu,j) = bbb(iu);
    end
% 	//---- set the source influence matrix for the wake sources
	for(i=1:p.n)
		for(j=p.n+1:p.n+p.nw) 
			p.dij(i,j) = p.bij(i,j);
        end
    end
	
% 	//**** now we need to calculate the influence of sources on the wake velocities
	
% 	//---- calculate dqtan/dgam and dqtan/dsig at the wake points

	for (i=p.n+1:p.n+p.nw)
		iw = i-p.n;
% 		//------ airfoil contribution at wake panel node
		[psi,psi_n,p]=psilin(i,p.x(i),p.y(i),p.nx(i),p.ny(i),false,true,p);
		for(j=1:p.n)
			p.cij(iw,j) = p.dqdg(j);
        end
		for(j=1:p.n) 
			p.dij(i,j) = p.dqdm(j);
        end
% 		//------ wake contribution
		[psi,psi_n,p]=pswlin(i,p.x(i),p.y(i),p.nx(i),p.ny(i),p);
		for(j=p.n+1: p.n+p.nw) 
			p.dij(i,j) = p.dqdm(j);
        end
    end
	
% 	//---- add on effect of all sources on airfoil vorticity which effects wake qtan
	for(i=p.n+1:p.n+p.nw)
		iw = i-p.n;
		
% 		//------ airfoil surface source contribution first
		for(j=1:p.n)
			sum = 0.0;
			for (k=1:p.n)
                sum = sum + p.cij(iw,k)*p.dij(k,j);
            end
			p.dij(i,j) = p.dij(i,j) + sum;
        end
		
% 		//------ wake source contribution next
		for(j=p.n+1:p.n+p.nw)
			sum = 0.0;
			for(k=1:p.n)
                sum = sum + p.cij(iw,k)*p.bij(k,j);
            end
			p.dij(i,j) = p.dij(i,j) + sum;
        end
		
    end
	
% 	//---- make sure first wake point has same velocity as trailing edge
	for(j=1:p.n+p.nw)
		p.dij(p.n+1,j) = p.dij(p.n,j);
    end
	
	
	p.lwdij = true;
% 	return true;
end
function [psi,psi_ni,p]=pswlin( i,  xi,  yi,  nxi,  nyi,p)
% ////--------------------------------------
% //     pswlin function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% 	//--------------------------------------------------------------------
% 	//	   calculates current streamfunction psi and tangential velocity
% 	//	   qtan at panel node or wake node i due to freestream and wake
% 	//	   sources p.sig.  also calculates sensitivity vectors dpsi/dsig
% 	//	   (p.dzdm) and dqtan/dsig (p.dqdm).
% 	//
% 	//			airfoil:  1   < i < p.n
% 	//			wake:	  p.n+1 < i < p.n+p.nw
% 	//--------------------------------------------------------------------
% 	
% 	 g1,g2,t1,t2;
% 	 x1i, x2i, yyi, x0,rs0,g0,t0;
% 	 dso, dsio, apan, rx1, rx2, ry1, ry2;
% 	 sx, sy, x1, x2, yy, rs1, rs2, sgn;
% 	 dxinv, psum, pdif, psx1, psx0, psyy, pdx1, pdx0, pdyy;
% 	 dsm, dsim, ssum, sdif, psni, pdni, psx2, pdx2, dsp, dsip;
% 	// nxi, nyi;
% 	int io,jo;
	io = i;
	
% 	p.cosa = cos(alfa);
% 	p.sina = sin(alfa);

	
	for(jo=p.n+1:p.n+p.nw)
		p.dzdm(jo) = 0.0;
		p.dqdm(jo) = 0.0;
    end
	
	psi	 = 0.0;
	psi_ni = 0.0;
	
	
	for(jo=p.n+1:p.n+p.nw-1)
		
		 jp = jo+1;
		 jm = jo-1;
		 jq = jp+1;
		if(jo==p.n+1) 
			jm = jo;
		
        else
			if(jo==p.n+p.nw-1) 
                jq = jp;
            end
        end
		dso = sqrt((p.x(jo)-p.x(jp))*(p.x(jo)-p.x(jp)) + (p.y(jo)-p.y(jp))* (p.y(jo)-p.y(jp)));
		dsio = 1.0 / dso;
		
		apan = p.apanel(jo);
		
		rx1 = xi - p.x(jo);
		ry1 = yi - p.y(jo);
		rx2 = xi - p.x(jp);
		ry2 = yi - p.y(jp);

		sx = (p.x(jp) - p.x(jo)) * dsio;
		sy = (p.y(jp) - p.y(jo)) * dsio;
		
		x1 = sx*rx1 + sy*ry1;
		x2 = sx*rx2 + sy*ry2;
		yy = sx*ry1 - sy*rx1;
		rs1 = rx1*rx1 + ry1*ry1;
		rs2 = rx2*rx2 + ry2*ry2;
		
		sgn =1.0;
		
		if(io>=p.n+1 && io<=p.n+p.nw) 
			sgn = 1.0;
        
        else
			sgn = 1*sign(yy);
        end
		
		if(io~=jo && rs1>0.0) 
			g1 = log(rs1);
			t1 = atan2(sgn*x1,sgn*yy) - (0.5- 0.5*sgn)*pi;
        
        else
			g1 = 0.0;
			t1 = 0.0;
        end
		
		if(io~=jp && rs2>0.0) 
			g2 = log(rs2);
			t2 = atan2(sgn*x2,sgn*yy) - (0.5- 0.5*sgn)*pi;
		
        else
			g2 = 0.0;
			t2 = 0.0;
        end
		x1i = sx*nxi + sy*nyi;
		x2i = sx*nxi + sy*nyi;
		yyi = sx*nyi - sy*nxi;
% 		//------- set up midpoint quantities
		x0 = 0.5*(x1+x2);
		rs0 = x0*x0 + yy*yy;
		g0 = log(rs0);
		t0 = atan2(sgn*x0,sgn*yy) - (0.5- 0.5*sgn)*pi;
		
% 		//------- calculate source contribution to psi	for  1-0  half-panel
		dxinv = 1.0/(x1-x0);
		psum = x0*(t0-apan) - x1*(t1-apan) + 0.5*yy*(g1-g0);
		pdif = ((x1+x0)*psum + rs1*(t1-apan) - rs0*(t0-apan)...
			+ (x0-x1)*yy) * dxinv;
		
		psx1 =  -(t1-apan);
		psx0 =    t0-apan;
		psyy =  0.5*(g1-g0);
		
		pdx1 = ((x1+x0)*psx1 + psum + 2.0*x1*(t1-apan) - pdif) * dxinv;
		pdx0 = ((x1+x0)*psx0 + psum - 2.0*x0*(t0-apan) + pdif) * dxinv;
		pdyy = ((x1+x0)*psyy + 2.0*(x0-x1 + yy*(t1-t0))	  ) * dxinv;
		
		dsm = sqrt((p.x(jp)-p.x(jm))*(p.x(jp)-p.x(jm)) + (p.y(jp)-p.y(jm))*(p.y(jp)-p.y(jm)));
		dsim = 1.0/dsm;
		
% 		////ccc		  sig0 = (p.sig(jp) - p.sig(jo))*dsio
% 		////ccc		  sig1 = (p.sig(jp) - p.sig(jm))*dsim
% 		////ccc		  ssum = sig0 + sig1
% 		////ccc		  sdif = sig0 - sig1

		
		ssum = (p.sig(jp) - p.sig(jo))*dsio + (p.sig(jp) - p.sig(jm))*dsim;
		sdif = (p.sig(jp) - p.sig(jo))*dsio - (p.sig(jp) - p.sig(jm))*dsim;
		
		psi = psi + p.qopi*(psum*ssum + pdif*sdif);
		
% 		//------- dpsi/dm
		p.dzdm(jm) = p.dzdm(jm) + p.qopi*(-psum*dsim + pdif*dsim);
		p.dzdm(jo) = p.dzdm(jo) + p.qopi*(-psum*dsio - pdif*dsio);
		p.dzdm(jp) = p.dzdm(jp) + p.qopi*( psum*(dsio+dsim) + pdif*(dsio-dsim));
		
% 		//------- dpsi/dni
		psni = psx1*x1i + psx0*(x1i+x2i)*0.5+ psyy*yyi;
		pdni = pdx1*x1i + pdx0*(x1i+x2i)*0.5+ pdyy*yyi;
		psi_ni = psi_ni + p.qopi*(psni*ssum + pdni*sdif);
		
		p.dqdm(jm) = p.dqdm(jm) + p.qopi*(-psni*dsim + pdni*dsim);
		p.dqdm(jo) = p.dqdm(jo) + p.qopi*(-psni*dsio - pdni*dsio);
		p.dqdm(jp) = p.dqdm(jp) + p.qopi*( psni*(dsio+dsim)+ pdni*(dsio-dsim));
		
		
% 		//------- calculate source contribution to psi	for  0-2  half-panel
		dxinv = 1.0/(x0-x2);
		psum = x2*(t2-apan) - x0*(t0-apan) + 0.5*yy*(g0-g2);
		pdif = ((x0+x2)*psum + rs0*(t0-apan) - rs2*(t2-apan)+ (x2-x0)*yy) * dxinv;
		
		psx0 =  -(t0-apan);
		psx2 =    t2-apan;
		psyy =  0.5*(g0-g2);
		
		pdx0 = ((x0+x2)*psx0 + psum + 2.0*x0*(t0-apan) - pdif) * dxinv;
		pdx2 = ((x0+x2)*psx2 + psum - 2.0*x2*(t2-apan) + pdif) * dxinv;
		pdyy = ((x0+x2)*psyy + 2.0*(x2-x0 + yy*(t0-t2))	  ) * dxinv;
		
		dsp = sqrt((p.x(jq)-p.x(jo))*(p.x(jq)-p.x(jo)) + (p.y(jq)-p.y(jo))*(p.y(jq)-p.y(jo)));
		dsip = 1.0/dsp;
		
% 		////ccc		  sig2 = (p.sig(jq) - p.sig(jo))*dsip
% 		////ccc		  sig0 = (p.sig(jp) - p.sig(jo))*dsio
% 		////ccc		  ssum = sig2 + sig0
% 		////ccc		  sdif = sig2 - sig0
		
		ssum = (p.sig(jq) - p.sig(jo))*dsip + (p.sig(jp) - p.sig(jo))*dsio;
		sdif = (p.sig(jq) - p.sig(jo))*dsip - (p.sig(jp) - p.sig(jo))*dsio;
		
		psi = psi + p.qopi*(psum*ssum + pdif*sdif);
		
% 		//------- dpsi/dm
		p.dzdm(jo) = p.dzdm(jo) + p.qopi*(-psum*(dsip+dsio)- pdif*(dsip-dsio));
		p.dzdm(jp) = p.dzdm(jp) + p.qopi*( psum*dsio - pdif*dsio);
		p.dzdm(jq) = p.dzdm(jq) + p.qopi*( psum*dsip + pdif*dsip);
		
% 		//------- dpsi/dni
		psni = psx0*(x1i+x2i)*0.5+ psx2*x2i + psyy*yyi;
		pdni = pdx0*(x1i+x2i)*0.5+ pdx2*x2i + pdyy*yyi;
		psi_ni = psi_ni + p.qopi*(psni*ssum + pdni*sdif);
		
		p.dqdm(jo) = p.dqdm(jo) + p.qopi*(-psni*(dsip+dsio)- pdni*(dsip-dsio));
		p.dqdm(jp) = p.dqdm(jp) + p.qopi*( psni*dsio - pdni*dsio);
		p.dqdm(jq) = p.dqdm(jq) + p.qopi*( psni*dsip + pdni*dsip);
		
	
    end
% 	return true;	  
end
function p=setbl(p)
% ////--------------------------------------
% //     setbl function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% %//-------------------------------------------------
% %//	   sets up the bl newton system coefficients for the current bl variables
% %//     and the edge velocities received from setup. the local bl system
% %//     coefficients are then incorporated into the global newton system.	
% %//-------------------------------------------------
% 
% 	int i, ibl,iv,iw, j, js,jv, jbl, is;
% 	int p.ile1,p.ile2,p.ite1,p.ite2,p.jvte1,p.jvte2;
% 	 p.usav(IVX+1,ISX);
% 	 p.u1_m(2*IVX+1), p.u2_m(2*IVX+1);
% 	 p.d1_m(2*IVX+1), p.d2_m(2*IVX+1);
% 	 p.ule1_m(2*IVX+1), p.ule2_m(2*IVX+1);
% 	 p.ute1_m(2*IVX+1), p.ute2_m(2*IVX+1);
% 	 p.msq_clmr, p.mdi;
% 	 p.herat,p.herat_ms;
% 	
% 	 p.clmr,ma_clmr,p.re_clmr;
% 	 p.ule1_a, p.ule2_a, p.u1_a, p.u2_a, p.d1_a, p.due1, p.due2, p.dds1, p.dds2;
% 	 xsi, cti, uei, thi, dsi, dswaki;
% 	 p.d2_a, p.d2_m2, p.d2_u2, p.dte_mte1, p.dte_ute1, p.dte_mte2, p.dte_ute2;
% 	 p.tte, p.cte, p.dte, p.dule1,p.dule2;
% 	 p.str, p.chx, p.chy, p.xtr, p.ytr, p.chsq;
% 	 p.xi_ule1, p.xi_ule2;

	 ami = 0.0;%//added arcds
	 p.tte_tte1= 0.0;%//added arcds
	 p.tte_tte2= 0.0;%//added arcds
	 p.cte_tte1= 0.0;%//added arcds
	 p.cte_tte2= 0.0;%//added arcds
	 p.cte_cte1= 0.0;%//added arcds
	 p.cte_cte2= 0.0;%//added arcds

	%//---- set the p.cl used to define mach, reynolds numbers
	if(p.lalfa) 
        p.clmr = p.cl;
    else
        p.clmr = p.clspec;
    end
	cti = 0.0; %//arcds added, otherwise variable is not initialized
	
	%//---- set current p.minf(p.cl)
	[p,ma_clmr, p.re_clmr]=mrcl(p.clmr,p);
	
	p.msq_clmr = 2.0*p.minf*ma_clmr;
	
	%//---- set compressibility parameter p.tklam and derivative tk_msq
	p=comset(p);
	
	%//---- set gas constant (= cp/cv)
	p.gambl = p.gamma;
	p.gm1bl = p.gamm1;
	
	%//---- set parameters for compressibility correction
	p.qinfbl  = p.qinf;
	p.tkbl	= p.tklam;
	p.tkbl_ms = p.tkl_msq;
	
	%//---- stagnation density and 1/enthalpy
	p.rstbl	 = ((1.0 + 0.5*p.gm1bl*p.minf*p.minf) ^(1.0/p.gm1bl));
	p.rstbl_ms = 0.5*p.rstbl/(1.0 + 0.5*p.gm1bl*p.minf*p.minf);
	p.hstinv	= p.gm1bl*(p.minf/p.qinfbl)*(p.minf/p.qinfbl) / (1.0 + 0.5*p.gm1bl*p.minf*p.minf);
	p.hstinv_ms = p.gm1bl*( 1.0/p.qinfbl)*( 1.0/p.qinfbl) / (1.0 + 0.5*p.gm1bl*p.minf*p.minf)...
		- 0.5*p.gm1bl*p.hstinv / (1.0 + 0.5*p.gm1bl*p.minf*p.minf);
	
	%//---- sutherland'p.s const./to	(assumes stagnation conditions are at stp)
	p.hvrat = 0.35;
	
	%//---- set reynolds number based on freestream density, velocity, viscosity
	p.herat	 = 1.0 - 0.5*p.qinfbl*p.qinfbl*p.hstinv;
	p.herat_ms =	 - 0.5*p.qinfbl*p.qinfbl*p.hstinv_ms;
	
	p.reybl	 = p.reinf * sqrt(p.herat*p.herat*p.herat) * (1.0+p.hvrat)/(p.herat+p.hvrat);
	p.reybl_re =		   sqrt(p.herat*p.herat*p.herat) * (1.0+p.hvrat)/(p.herat+p.hvrat);
	p.reybl_ms = p.reybl * (1.5/p.herat - 1.0/(p.herat+p.hvrat))*p.herat_ms;
	
	p.amcrit = p.acrit;
	
	%//---- save te thickness
	p.dwte = p.wgap(1);
	
	if(p.lblini==false)
		%//----- initialize bl by marching with ue (fudge at separation)
		disp(' initializing bl ...');
% 		CString p.str;
% 		p.str.Format(" Initializing bl ...\r\p.n");
% 		if(m_bTrace)pXFile->WriteString(p.str);
        
		p=mrchue(p);
		p.lblini = true;
    end

	%//---- march bl with current ue and ds to establish transition
	p=mrchdu(p);
	
	for (is=1:2)
		for(ibl=2:p.nbl(is)) 
			p.usav(ibl,is) = p.uedg(ibl,is);
        end
    end
	
	p=ueset(p);
	
	for (is=1:2)
		for(ibl=2:p.nbl(is))
			 temp = p.usav(ibl,is);
			p.usav(ibl,is) = p.uedg(ibl,is);
			p.uedg(ibl,is) = temp;
        end
    end
	p.ile1 = p.ipan(2,1);
	p.ile2 = p.ipan(2,2);
	p.ite1 = p.ipan(p.iblte(1),1);
	p.ite2 = p.ipan(p.iblte(2),2);
	
	p.jvte1 = p.isys(p.iblte(1),1);
	p.jvte2 = p.isys(p.iblte(2),2);
	
	p.dule1 = p.uedg(2,1) - p.usav(2,1);
	p.dule2 = p.uedg(2,2) - p.usav(2,2);
	
	%//---- set le and te ue sensitivities wrt all m values
	for(js=1:2)
		for(jbl=2:p.nbl(js))
			j  = p.ipan(jbl,js);
			jv = p.isys(jbl,js);
			p.ule1_m(jv) = -p.vti(	     2,1)*p.vti(jbl,js)*p.dij(p.ile1,j);
			p.ule2_m(jv) = -p.vti(	     2,2)*p.vti(jbl,js)*p.dij(p.ile2,j);
			p.ute1_m(jv) = -p.vti(p.iblte(1),1)*p.vti(jbl,js)*p.dij(p.ite1,j);
			p.ute2_m(jv) = -p.vti(p.iblte(2),2)*p.vti(jbl,js)*p.dij(p.ite2,j);
        end
    end
	
	p.ule1_a = p.uinv_a(2,1);
	p.ule2_a = p.uinv_a(2,2);
	
% 	CString str1;
% 	str1.Format(" \r\p.n");
% 	if(m_bTrace)pXFile->WriteString(str1);
	%//*** go over each boundary layer/p.wake
	for(is=1:2)			
		%//---- there is no station "1" at similarity, so zero everything out
		for(js=1:2)
			for(jbl=2:p.nbl(js))
				jv = p.isys(jbl,js);
				p.u1_m(jv) = 0.0;
				p.d1_m(jv) = 0.0;
            end
        end
		p.u1_a = 0.0;
		p.d1_a = 0.0;
		
		p.due1 = 0.0;
		p.dds1 = 0.0;
		
		%//---- similarity station pressure gradient parameter  x/u du/dx
		ibl = 2;
		p.bule = 1.0;
		
		%//---- set forced transition arc length position
		p=xifset(is,p);
		
		p.tran = false;
		p.turb = false;
		
		%//**** sweep downstream setting up bl equation linearizations
		for(ibl=2:p.nbl(is))
			
			iv	= p.isys(ibl,is);
			
			p.simi = (ibl==2);
			p.wake = (ibl>p.iblte(is));
			p.tran = (ibl==p.itran(is));
			p.turb = (ibl>p.itran(is));
			
			i = p.ipan(ibl,is);

			%//---- set primary variables for current station
			xsi = p.xssi(ibl,is);
			if(ibl<p.itran(is)) 
                ami = p.ctau(ibl,is);
            else
                cti = p.ctau(ibl,is);
            end
			uei = p.uedg(ibl,is);
			thi = p.thet(ibl,is);
			p.mdi = p.mass(ibl,is);

			dsi = p.mdi/uei;
			
			if(p.wake) 
				iw = ibl - p.iblte(is);
				dswaki = p.wgap(iw);
			
            else
                dswaki = 0.0;
            end
			
			%//---- set derivatives of dsi (= d2)
			p.d2_m2 =  1.0/uei;
			p.d2_u2 = -dsi /uei;
			
			for(js=1:2)
				for(jbl=2:p.nbl(js))
					j  = p.ipan(jbl,js);
					jv = p.isys(jbl,js);
					p.u2_m(jv) = -p.vti(ibl,is)*p.vti(jbl,js)*p.dij(i,j);
					p.d2_m(jv) = p.d2_u2*p.u2_m(jv);
                end
            end
			p.d2_m(iv) = p.d2_m(iv) + p.d2_m2;
			
			p.u2_a = p.uinv_a(ibl,is);
			p.d2_a = p.d2_u2*p.u2_a;
			
			%//---- "forced" changes due to mismatch between p.uedg and p.usav=uinv+p.dij*p.mass
			p.due2 = p.uedg(ibl,is) - p.usav(ibl,is);
			p.dds2 = p.d2_u2*p.due2;
			
			p=blprv(xsi,ami,cti,thi,dsi,dswaki,uei,p);
			p=blkin(p);
			
			%//---- check for transition and set p.tran, p.xt, etc. if found
			if(p.tran) 
				p=trchek(p);
				ami = p.ampl2;
            end

			if(ibl==p.itran(is) && p.tran==false)
				disp('setbl: p.xtr???  ');
% 				CString p.str;
% 				p.str.Format("setbl: p.xtr???  n1=%d n2=%d: \r\p.n", ampl1, ampl2);
% 				if(m_bTrace)pXFile->WriteString(p.str);

            end
			
			%//---- assemble 10x4 linearized system for dctau, dth, dds, due, dxi
			%//	   at the previous "1" station and the current "2" station
			
			if(ibl==p.iblte(is)+1) 
				
				%//----- define quantities at start of p.wake, adding te base thickness to dstar
				p.tte = p.thet(p.iblte(1),1) + p.thet(p.iblte(2),2);
				p.dte = p.dstr(p.iblte(1),1) + p.dstr(p.iblte(2),2) + p.ante;
				p.cte = ( p.ctau(p.iblte(1),1)*p.thet(p.iblte(1),1) + p.ctau(p.iblte(2),2)*p.thet(p.iblte(2),2) ) / p.tte;
				p=tesys(p.cte,p.tte,p.dte,p);
				
				p.tte_tte1 = 1.0;
				p.tte_tte2 = 1.0;
				p.dte_mte1 =				1.0  / p.uedg(p.iblte(1),1);
				p.dte_ute1 = -p.dstr(p.iblte(1),1) / p.uedg(p.iblte(1),1);
				p.dte_mte2 =				1.0  / p.uedg(p.iblte(2),2);
				p.dte_ute2 = -p.dstr(p.iblte(2),2) / p.uedg(p.iblte(2),2);
				p.cte_cte1 = p.thet(p.iblte(1),1)/p.tte;
				p.cte_cte2 = p.thet(p.iblte(2),2)/p.tte;
				p.cte_tte1 = (p.ctau(p.iblte(1),1) - p.cte)/p.tte;
				p.cte_tte2 = (p.ctau(p.iblte(2),2) - p.cte)/p.tte;
				
				%//----- re-define d1 sensitivities wrt m since d1 depends on both te ds values
				for (js=1:2)
					for (jbl=2:p.nbl(js))
						j  = p.ipan(jbl,js);
						jv = p.isys(jbl,js);
						p.d1_m(jv) = p.dte_ute1*p.ute1_m(jv) + p.dte_ute2*p.ute2_m(jv);
                    end
                end
				p.d1_m(p.jvte1) = p.d1_m(p.jvte1) + p.dte_mte1;
				p.d1_m(p.jvte2) = p.d1_m(p.jvte2) + p.dte_mte2;
				
				%//----- "forced" changes from  p.uedg --- p.usav=uinv+p.dij*p.mass	mismatch
				p.due1 = 0.0;
				p.dds1 = p.dte_ute1*(p.uedg(p.iblte(1),1) - p.usav(p.iblte(1),1))...
					+ p.dte_ute2*(p.uedg(p.iblte(2),2) - p.usav(p.iblte(2),2));
			
            else
				p=blsys(p);
            end
			
			
			%//---- save wall shear and equil. max shear coefficient for plotting output
			p.tau(ibl,is) = 0.5*p.r2*p.u2*p.u2*p.cf2;
			p.dis(ibl,is) =	 p.r2*p.u2*p.u2*p.u2*p.di2*p.hs2*0.5;
			p.ctq(ibl,is) = p.cq2;
			p.delt(ibl,is) = p.de2;
			p.uslp(ibl,is) = 1.60/(1.0+p.us2);
			
			%//---- set xi sensitivities wrt le ue changes
			if(is==1) 
				p.xi_ule1 =  p.sst_go;
				p.xi_ule2 = -p.sst_gp;
			
            else
				p.xi_ule1 = -p.sst_go;
				p.xi_ule2 =  p.sst_gp;
            end
			
			%//---- stuff bl system coefficients into main jacobian matrix
			
			for( jv=1:p.nsys)
				p.vm(1,jv,iv) = p.vs1(1,3)*p.d1_m(jv) + p.vs1(1,4)*p.u1_m(jv)...
            		+ p.vs2(1,3)*p.d2_m(jv) + p.vs2(1,4)*p.u2_m(jv)...
					+ (p.vs1(1,5) + p.vs2(1,5) + p.vsx(1))...
					*(p.xi_ule1*p.ule1_m(jv) + p.xi_ule2*p.ule2_m(jv));
            end
			
			p.vb(1,1,iv) = p.vs1(1,1);
			p.vb(1,2,iv) = p.vs1(1,2);
			
			p.va(1,1,iv) = p.vs2(1,1);
			p.va(1,2,iv) = p.vs2(1,2);
			
			if(p.lalfa)  
                p.vdel(1,2,iv) = p.vsr(1)*p.re_clmr + p.vsm(1)*p.msq_clmr;
            else
                p.vdel(1,2,iv) =(p.vs1(1,4)*p.u1_a + p.vs1(1,3)*p.d1_a)...
				+ (p.vs2(1,4)*p.u2_a + p.vs2(1,3)*p.d2_a)...
				+ (p.vs1(1,5) + p.vs2(1,5) + p.vsx(1))...
				*(p.xi_ule1*p.ule1_a + p.xi_ule2*p.ule2_a);
            end
			
			p.vdel(1,1,iv) = p.vsrez(1)...
				+ (p.vs1(1,4)*p.due1 + p.vs1(1,3)*p.dds1)...
				+ (p.vs2(1,4)*p.due2 + p.vs2(1,3)*p.dds2)...
				+ (p.vs1(1,5) + p.vs2(1,5) + p.vsx(1))...
				*(p.xi_ule1*p.dule1 + p.xi_ule2*p.dule2);
			
			for(jv=1:p.nsys)
				p.vm(2,jv,iv) = p.vs1(2,3)*p.d1_m(jv) + p.vs1(2,4)*p.u1_m(jv)...
					+ p.vs2(2,3)*p.d2_m(jv) + p.vs2(2,4)*p.u2_m(jv)...
					+ (p.vs1(2,5) + p.vs2(2,5) + p.vsx(2))...
					*(p.xi_ule1*p.ule1_m(jv) + p.xi_ule2*p.ule2_m(jv));
            end
			p.vb(2,1,iv)	= p.vs1(2,1);
			p.vb(2,2,iv)	= p.vs1(2,2);
			
			p.va(2,1,iv) = p.vs2(2,1);
			p.va(2,2,iv) = p.vs2(2,2);
			
			if(p.lalfa) p.vdel(2,2,iv) = p.vsr(2)*p.re_clmr + p.vsm(2)*p.msq_clmr;
			else		 p.vdel(2,2,iv) = ...
				(p.vs1(2,4)*p.u1_a + p.vs1(2,3)*p.d1_a)...
				+ (p.vs2(2,4)*p.u2_a + p.vs2(2,3)*p.d2_a)...
				+ (p.vs1(2,5) + p.vs2(2,5) + p.vsx(2))...
				*(p.xi_ule1*p.ule1_a + p.xi_ule2*p.ule2_a);
            end
			
			p.vdel(2,1,iv) = p.vsrez(2)...
				+ (p.vs1(2,4)*p.due1 + p.vs1(2,3)*p.dds1)...
				+ (p.vs2(2,4)*p.due2 + p.vs2(2,3)*p.dds2)...
				+ (p.vs1(2,5) + p.vs2(2,5) + p.vsx(2))...
				*(p.xi_ule1*p.dule1 + p.xi_ule2*p.dule2);
			

% 			//memory overlap problem
			for(jv=1:p.nsys)
				p.vm(3,jv,iv) = p.vs1(3,3)*p.d1_m(jv) + p.vs1(3,4)*p.u1_m(jv)...
							  + p.vs2(3,3)*p.d2_m(jv) + p.vs2(3,4)*p.u2_m(jv)...
							  + (p.vs1(3,5) + p.vs2(3,5) + p.vsx(3))...
								*(p.xi_ule1*p.ule1_m(jv) + p.xi_ule2*p.ule2_m(jv));
            end
			
			p.vb(3,1,iv) = p.vs1(3,1);
			p.vb(3,2,iv) = p.vs1(3,2);
			
			p.va(3,1,iv) = p.vs2(3,1);
			p.va(3,2,iv) = p.vs2(3,2);
			
			if(p.lalfa) 
                p.vdel(3,2,iv) = p.vsr(3)*p.re_clmr + p.vsm(3)*p.msq_clmr;
			else		 p.vdel(3,2,iv) = ...
				(p.vs1(3,4)*p.u1_a + p.vs1(3,3)*p.d1_a)...
				+ (p.vs2(3,4)*p.u2_a + p.vs2(3,3)*p.d2_a)...
				+ (p.vs1(3,5) + p.vs2(3,5) + p.vsx(3))...
				*(p.xi_ule1*p.ule1_a + p.xi_ule2*p.ule2_a);
            end
			
			p.vdel(3,1,iv) = p.vsrez(3)...
				+ (p.vs1(3,4)*p.due1 + p.vs1(3,3)*p.dds1)...
				+ (p.vs2(3,4)*p.due2 + p.vs2(3,3)*p.dds2)...
				+ (p.vs1(3,5) + p.vs2(3,5) + p.vsx(3))...
				*(p.xi_ule1*p.dule1 + p.xi_ule2*p.dule2);
			
			if(ibl==p.iblte(is)+1) 
				
				%//----- redefine coefficients for p.tte, p.dte, etc
				p.vz(1,1)	= p.vs1(1,1)*p.cte_cte1;
				p.vz(1,2)	= p.vs1(1,1)*p.cte_tte1 + p.vs1(1,2)*p.tte_tte1;
				p.vb(1,1,iv) = p.vs1(1,1)*p.cte_cte2;
				p.vb(1,2,iv) = p.vs1(1,1)*p.cte_tte2 + p.vs1(1,2)*p.tte_tte2;
				
				p.vz(2,1)	= p.vs1(2,1)*p.cte_cte1;
				p.vz(2,2)	= p.vs1(2,1)*p.cte_tte1 + p.vs1(2,2)*p.tte_tte1;
				p.vb(2,1,iv) = p.vs1(2,1)*p.cte_cte2;
				p.vb(2,2,iv) = p.vs1(2,1)*p.cte_tte2 + p.vs1(2,2)*p.tte_tte2;
				
				p.vz(3,1)	= p.vs1(3,1)*p.cte_cte1;
				p.vz(3,2)	= p.vs1(3,1)*p.cte_tte1 + p.vs1(3,2)*p.tte_tte1;
				p.vb(3,1,iv) = p.vs1(3,1)*p.cte_cte2;
				p.vb(3,2,iv) = p.vs1(3,1)*p.cte_tte2 + p.vs1(3,2)*p.tte_tte2;
				
            end
			
			%//---- turbulent intervals will follow if currently at transition interval
			if(p.tran) 
				p.turb = true;
				
				%//------ save transition location
				p.itran(is) = ibl;
				tforce(is) = p.trforc;
				p.xssitr(is) = p.xt;
				
				%//------ interpolate airfoil geometry to find transition x/c
				%//		(for user output)
				if(is==1) 
                    p.str = p.sst - p.xt;
                else
                    p.str = p.sst + p.xt;
                end
				p.chx = p.xte - p.xle;
				p.chy = p.yte - p.yle;
				p.chsq = p.chx*p.chx + p.chy*p.chy;
				p.xtr = seval(p.str,p.x,p.xp,p.s,p.n);
				p.ytr = seval(p.str,p.y,p.yp,p.s,p.n);
				p.xoctr(is) = ((p.xtr-p.xle)*p.chx + (p.ytr-p.yle)*p.chy)/p.chsq;
				p.yoctr(is) = ((p.ytr-p.yle)*p.chx - (p.xtr-p.xle)*p.chy)/p.chsq;
            end
			
			p.tran = false;
			
			if(ibl==p.iblte(is)) 
				%//----- set "2" variables at te to p.wake correlations for next station
				
				p.turb = true;
				p.wake = true;
				p=blvar(3,p);
				p=blmid(3,p);
            end
			
			for(js=1:2)
				for(jbl=2:p.nbl(js))
					jv = p.isys(jbl,js);
					p.u1_m(jv) = p.u2_m(jv);
					p.d1_m(jv) = p.d2_m(jv);
                end
            end
			
			p.u1_a = p.u2_a;
			p.d1_a = p.d2_a;
			
			p.due1 = p.due2;
			p.dds1 = p.dds2;
			
			%//---- set bl variables for next station
%//			for (int icom=1; icom<= ncom;icom++)	com1(icom) = com2(icom);
			p=stepbl(p);
			
			%//---- next streamwise station
        end
			
% 		CString strOut;
		if(p.tforce(is))
			disp('     Side %d, forced transition ');
% 			//TRACE(strOut);
% 			if(m_bTrace)pXFile->WriteString(strOut);

		
        else
			disp('     Side %d,  free  transition ');
% 			//TRACE(strOut);
% 			if(m_bTrace)pXFile->WriteString(strOut);
%             end
        end
		
		%//---- next airfoil side
    end
% 	return true;
end
function p=mrchue(p)
% ////--------------------------------------
% //     mrchue function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% 	%%//----------------------------------------------------
% 	%%//     marches the bls and p.wake in direct mode using
% 	%%//     the p.uedg array. if separation is encountered,
% 	%%//     a plausible value of hk extrapolated from
% 	%%//     upstream is prescribed instead.  continuous
% 	%%//     checking of transition onset is performed.
% 	%%//----------------------------------------------------
% 	
% 	CString str;
% 	bool direct;
% 	int ibl, ibm, iw;
% 	 msq, ratlen,dsw,hklim;
% 	 hlmax, htmax, xsi, uei, ucon, tsq, thi, ami, cti, dsi;
% 	 dswaki;
% 	  htest, hktest, dummy;
% 	 cst;			%//added arcds
	 cte = 0.0;	%//added arcds
	 dte = 0.0;	%//added arcds
	 tte = 0.0;	%//added arcds
	 dmax = 0.0;	%//added arcds
	 hmax = 0.0;	%//added arcds
	 p.htarg = 0.0;	%//added arcds

    goto=0;
	%%//---- shape parameters for separation criteria
	hlmax = 3.8;
	htmax = 2.5;
	
	for (is=1:2)%//2000
		
% 		CString str;
		disp('    Side %d ...\r\n');
% 		if(m_bTrace)pXFile->WriteString(str);

		%%//---- set forced transition arc length position
		p=xifset(is,p);
		
		%%//---- initialize similarity station with thwaites' formula
		%%//	ibl = 2;
		xsi = p.xssi(2,is);
		uei = p.uedg(2,is);


		%%//      bule = log(p.uedg(ibl+1,is)/uei) / log(p.xssi(ibl+1,is)/xsi)
		%%//      bule = max( -.08 , bule )
		bule = 1.0;
		ucon = uei/(xsi^bule);
		tsq = 0.45/(ucon*(5.0*bule+1.0)*p.reybl) * (xsi^(1.0-bule));
		thi = sqrt(tsq);
		dsi = 2.2*thi;
		ami = 0.0;
		
		%%//---- initialize p.ctau for first turbulent station
		cti = 0.03;
		
		p.tran = false;
		p.turb = false;
		p.itran(is) = p.iblte(is);
		
		%%//---- march downstream
		for (ibl=2:p.nbl(is))%%// 1000
			ibm = ibl-1;
			iw = ibl - p.iblte(is);
			p.simi = (ibl==2);
			p.wake = ibl>p.iblte(is);
			
			%%//------ prescribed quantities
			xsi = p.xssi(ibl,is);
			uei = p.uedg(ibl,is);

		
			if(p.wake) 
				iw = ibl - p.iblte(is);
				dswaki = p.wgap(iw);
			
			else dswaki = 0.0;
            end
			
			direct = true;
			
			%%//------ newton iteration loop for current station 
			for ( itbl=1:25)%//100
				
				%%//-------- assemble 10x3 linearized system for dctau, dth, dds, due, dxi
				%%//         at the previous "1" station and the current "2" station
				%%//         (the "1" station coefficients will be ignored)
				
				p=blprv(xsi,ami,cti,thi,dsi,dswaki,uei,p);
				p=blkin(p);
				
				%%//-------- check for transition and set appropriate flags and things
				if((p.simi==false) && (p.turb==false)) 
					p=trchek(p);
					ami = p.ampl2;
					
					%%//--------- fixed bug   md 7 jun 99
					if(p.tran) 
						p.itran(is) = ibl;
						if(cti<=0.0) 
							cti = 0.03;
							p.s2 = cti;
                        end
					
                    else
                        p.itran(is) = ibl+2;
                    end
                end
				
				if(ibl==p.iblte(is)+1) 
					tte = p.thet(p.iblte(1),1) + p.thet(p.iblte(2),2);
					dte = p.dstr(p.iblte(1),1) + p.dstr(p.iblte(2),2) + p.ante;
					cte = ( p.ctau(p.iblte(1),1)*p.thet(p.iblte(1),1)...
						+ p.ctau(p.iblte(2),2)*p.thet(p.iblte(2),2) ) / tte;
					p=tesys(cte,tte,dte,p);
				
				else 
					p=blsys(p);
                end
				if(direct) 
					%%//--------- try direct mode (set due = 0 in currently empty 4th line)
					p.vs2(4,1) = 0.0;
					p.vs2(4,2) = 0.0;
					p.vs2(4,3) = 0.0;
					p.vs2(4,4) = 1.0;
					p.vsrez(4) = 0.0;
					%%//--------- solve newton system for current "2" station
					p.vsrez=Gauss(4,p.vs2,p.vsrez);
					%%//--------- determine max changes and underrelax if necessary
					dmax = max( abs(p.vsrez(2)/thi), abs(p.vsrez(3)/dsi) );
					if(ibl<p.itran(is)) 
                        dmax = max(dmax,abs(p.vsrez(1)/10.0));
                    end
					if(ibl>=p.itran(is))
                        dmax = max(dmax,abs(p.vsrez(1)/cti ));
                    end
					rlx = 1.0;
					if(dmax>0.3) 
                        rlx = 0.3/dmax;
                    end
					%%//--------- see if direct mode is not applicable
					if(ibl ~= p.iblte(is)+1) 
						%%//---------- calculate resulting kinematic shape parameter hk
						msq = uei*uei*p.hstinv / (p.gm1bl*(1.0 - 0.5*uei*uei*p.hstinv));
						htest = (dsi + rlx*p.vsrez(3)) / (thi + rlx*p.vsrez(2));
						[hktest, dummy, dummy]=hkin(htest, msq);
						
						%%//---------- decide whether to do direct or inverse problem based on hk
						if(ibl<p.itran(is))
                            hmax = hlmax;
                        end
						if(ibl>=p.itran(is)) 
                            hmax = htmax;
                        end
						direct = (hktest<hmax);
                    end
					if(direct) 
						%%//---------- update as usual
						if(ibl>=p.itran(is))
                            cti = cti + rlx*p.vsrez(1);
                        end
						thi = thi + rlx*p.vsrez(2);
						dsi = dsi + rlx*p.vsrez(3);
					
                    else
						%%//---------- set prescribed hk for inverse calculation at the current station
						if(ibl<p.itran(is)) 
							%%//----------- laminar case: relatively slow increase in hk downstream
							p.htarg = p.hk1 + 0.03*(p.x2-p.x1)/p.t1;
                        elseif(ibl==p.itran(is)) 
							%%//----------- transition interval: weighted laminar and turbulent case
							p.htarg = p.hk1 + (0.03*(p.xt-p.x1) - 0.15*(p.x2-p.xt))/p.t1;
						
                        elseif(p.wake) 
							%%//----------- turbulent p.wake case:
							%%//--          asymptotic p.wake behavior with approximate backward euler
							cst = 0.03*(p.x2-p.x1)/p.t1;
							p.hk2 = p.hk1;
							p.hk2 = p.hk2 - (p.hk2 +     cst*(p.hk2-1.0)*(p.hk2-1.0)*(p.hk2-1.0) - p.hk1)...
								/(1.0 + 3.0*cst*(p.hk2-1.0)*(p.hk2-1.0));
							p.hk2 = p.hk2 - (p.hk2 +     cst*(p.hk2-1.0)*(p.hk2-1.0)*(p.hk2-1.0) - p.hk1)...
								/(1.0 + 3.0*cst*(p.hk2-1.0)*(p.hk2-1.0));
							p.hk2 = p.hk2 - (p.hk2 +     cst*(p.hk2-1.0)*(p.hk2-1.0)*(p.hk2-1.0) - p.hk1)...
								/(1.0 + 3.0*cst*(p.hk2-1.0)*(p.hk2-1.0));
							p.htarg = p.hk2;
						
                        else
                                p.htarg = p.hk1 - 0.15*(p.x2-p.x1)/p.t1;%%//----------- turbulent case: relatively fast decrease in hk downstream
                        end
						%%//---------- limit specified hk to something reasonable
						if(p.wake) 
                            p.htarg = max(p.htarg , 1.01);
						else p.htarg = max(p.htarg , hmax);
                        end
% 						CString str;
						disp('     mrchue: inverse mode at ');
% 						if(m_bTrace)pXFile->WriteString(str);


						%%//---------- try again with prescribed hk

						goto=100;
                        break
                    end
                
                else
					%%//-------- inverse mode (force hk to prescribed value p.htarg)
					p.vs2(4,1) = 0.0;
					p.vs2(4,2) = hk2_t2;
					p.vs2(4,3) = hk2_d2;
					p.vs2(4,4) = hk2_u2;
					p.vsrez(4) = p.htarg - p.hk2;
					p.vsrez=Gauss(4,p.vs2,p.vsrez);

					dmax = max( abs(p.vsrez(2)/thi),abs(p.vsrez(3)/dsi)  );
					if(ibl>=p.itran(is)) 
                        dmax = max( dmax , abs(p.vsrez(1)/cti));
                    end
					rlx = 1.0;
					if(dmax>0.3)
                        rlx = 0.3/dmax;
                    end
					%%//--------- update variables
					if(ibl>=p.itran(is)) 
                        cti = cti + rlx*p.vsrez(1);
                    end
					thi = thi + rlx*p.vsrez(2);
					dsi = dsi + rlx*p.vsrez(3);
					uei = uei + rlx*p.vsrez(4);

                end
				%%//-------- eliminate absurd transients
              if goto~=100
				if(ibl>=p.itran(is)) 
					cti = min(cti, 0.30);
					cti = max(cti, 0.0000001);
                end
				if(ibl<=p.iblte(is))  
                    hklim = 1.02;
                else
                    hklim = 1.00005;
                end
				msq = uei*uei*p.hstinv / (p.gm1bl*(1.0 - 0.5*uei*uei*p.hstinv));
				dsw = dsi - dswaki;
				dslim(dsw,thi,msq,hklim);
				dsi = dsw + dswaki;
				if(dmax<=0.00001) 
                    goto=110;
                    break
                end
              end
% stop100:
% 			int nothing;
			nothing = 1;
			
            end%//end itbl loop
           if goto~=110

			%//TRACE(" mrchue: convergence failed at %d,  side %d, res = %f\n", ibl, is, dmax);
% 			disp('     mrchue: convergence failed at %d,  side %d, res = %.3f\r\n');
% 			if(m_bTrace)pXFile->WriteString(str);

			%%//------ the current unconverged solution might still be reasonable...
			if(dmax > 0.1) 
				
				%%//------- the current solution is garbage --> extrapolate values instead
				if(ibl>3) 
					if(ibl<=p.iblte(is)) 
						thi = p.thet(ibm,is) * sqrt(p.xssi(ibl,is)/p.xssi(ibm,is));
						dsi = p.dstr(ibm,is) * sqrt(p.xssi(ibl,is)/p.xssi(ibm,is));
                    
                    else
						if(ibl==p.iblte(is)+1) 
							cti = cte;
							thi = tte;
							dsi = dte;
						
                        else
							thi = p.thet(ibm,is);
							ratlen = (p.xssi(ibl,is)-p.xssi(ibm,is)) / (10.0*p.dstr(ibm,is));
							dsi = (p.dstr(ibm,is) + thi*ratlen) / (1.0 + ratlen);
                        end
                    end
					if(ibl==p.itran(is))
                        cti = 0.05;
                    end
					if(ibl>p.itran(is)) 
                        cti = p.ctau(ibm,is);
                    end
					
					uei = p.uedg(ibl,is);

					if(ibl>2 && ibl<p.nbl(is)) 
                        uei = 0.5*(p.uedg(ibl-1,is) + p.uedg(ibl+1,is));
                    end
                end
            end
			%//109
			p=blprv(xsi,ami,cti,thi,dsi,dswaki,uei,p);
			p=blkin(p);
			%%//------- check for transition and set appropriate flags and things
			if((p.simi==false) && (p.turb==false)) 
				p=trchek(p);
				ami = p.ampl2;
				if(	  p.tran) 
                    p.itran(is) = ibl;
                end
				if(p.tran==false)
                    p.itran(is) = ibl+2;
                end
            end
			%%//------- set all other extrapolated values for current station
			if(ibl<p.itran(is)) 
                p=blvar(1,p);
            end
			if(ibl>=p.itran(is))
                p=blvar(2,p);
            end
			if(p.wake)  
                p=blvar(3,p);
			end
            if(ibl<p.itran(is))  
                p=blmid(1,p);
                end
			if(ibl>=p.itran(is)) 
                p=blmid(2,p);
                end
			if(p.wake)
                p=blmid(3,p);
                end
			%%//------ pick up here after the newton iterations
           end
			%%//------ store primary variables
			if(ibl<p.itran(is))
                p.ctau(ibl,is) = ami;
            end
			if(ibl>=p.itran(is))
                p.ctau(ibl,is) = cti;
            end
			p.thet(ibl,is) = thi;
			p.dstr(ibl,is) = dsi;
			p.uedg(ibl,is) = uei;
			p.mass(ibl,is) = dsi*uei;
			p.tau(ibl,is)  = 0.5*p.r2*p.u2*p.u2*p.cf2;
			p.dis(ibl,is)  = 	p.r2*p.u2*p.u2*p.u2*p.di2*p.hs2*0.5;
			p.ctq(ibl,is)  = p.cq2;
			p.delt(ibl,is) = p.de2;
			
			%%//------ set "1" variables to "2" variables for next streamwise station
			p=blprv(xsi,ami,cti,thi,dsi,dswaki,uei,p);
			p=blkin(p);

			p=stepbl(p);
			
			
			%%//------ turbulent intervals will follow transition interval or te
			if(p.tran || ibl==p.iblte(is)) 
				p.turb = true;
				
				%%//------- save transition location
				p.tforce(is) = p.trforc;
				p.xssitr(is) = p.xt;
            end
			
			p.tran = false;
			
% 			if(ibl==p.iblte(is)) 
% 				thi = p.thet(p.iblte(1),1) + p.thet(p.iblte(2),2);
% 				dsi = p.dstr(p.iblte(1),1) + p.dstr(p.iblte(2),2) + p.ante;
%             end
        end%%// 1000 continue : end ibl loop
    end%%// 2000 continue : end is loop
% 	return true;
end
function p=tesys( cte,  tte,  dte,p)
% ////--------------------------------------
% //     tesys function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% 	//--------------------------------------------------------
% 	//	   sets up "dummy" bl system between airfoil te point 
% 	//	   and first wake point infinitesimally behind te.
% 	//--------------------------------------------------------
	
	for( k=1: 4)
		p.vsrez(k) = 0.0;
		p.vsm(k)	 = 0.0;
		p.vsr(k)	 = 0.0;
		p.vsx(k)	 = 0.0;
		for (l=1:5)
			p.vs1(k,l) = 0.0;
			p.vs2(k,l) = 0.0;
        end
    end
	
	p=blvar(3,p);
	
	p.vs1(1,1) = -1.0;
	p.vs2(1,1) = 1.0;
	p.vsrez(1) = cte - p.s2;
	
	p.vs1(2,2) = -1.0;
	p.vs2(2,2) = 1.0;
	p.vsrez(2) = tte - p.t2;
	
	p.vs1(3,3) = -1.0;
	p.vs2(3,3) = 1.0;
	p.vsrez(3) = dte - p.d2 - p.dw2;
	
% 	return true;
% }
end
function [ax,  ax_hk,  ax_th,  ax_rt]=dampl( hk,  th,  rt)
% ////--------------------------------------
% //     dampl function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

%//==============================================================
%//     amplification rate routine for envelope e^n method.
%//     reference:
%//                drela, m., giles, m.,
%//               "viscous/inviscid analysis of transonic and
%//                low reynolds number airfoils",
%//                aiaa journal, oct. 1987.
%//
%//     new version.   march 1991       (latest bug fix  july 93)
%//          - m(h) correlation made more accurate up to h=20
%//          - for h > 5, non-similar profiles are used 
%//            instead of falkner-skan profiles.  these 
%//            non-similar profiles have smaller reverse 
%//            velocities, are more representative of typical 
%//            separation bubble profiles.
%//--------------------------------------------------------------
%//
%//     input :   hk     kinematic shape parameter
%//               th     momentum thickness
%//               rt     momentum-thickness reynolds number
%//
%//     output:   ax     envelope spatial amplification rate
%//               ax_(.) sensitivity of ax to parameter (.)
%//
%//
%//     usage: the log of the envelope amplitude n(x) is
%//            calculated by integrating ax (= dn/dx) with
%//            respect to the streamwise distance x.
%//                      x
%//                     /
%//              n(x) = | ax(h(x),th(x),rth(x)) dx
%//                     /
%//                      0
%//            the integration can be started from the leading
%//            edge since ax will be returned as zero when rt
%//            is below the critical rtheta.  transition occurs
%//            when n(x) reaches ncrit (ncrit= 9 is "standard").
%//==============================================================
  %//    implicit real (a-h,m,o-z)
	 dgr = 0.08;

% 	 hmi, hmi_hk,aa, aa_hk, bb, bb_hk, grcrit, grc_hk, gr, gr_rt;
% 	 rnorm, rn_hk, rn_rt, rfac, rfac_hk, rfac_rt;
% 	 rfac_rn, arg_hk,ex, f_arg,	ex_hk;
% 	 af, af_hmi, af_hk, dadr, dadr_hk;
	hmi = 1.0/(hk - 1.0);
	hmi_hk = -hmi*hmi;
	
	%//---- log10(critical rth) - h   correlation for falkner-skan profiles
	aa    = 2.492*(hmi^0.43);
	aa_hk =   (aa/hmi)*0.43 * hmi_hk;
	bb    = tanh(14.0*hmi - 9.24);
	bb_hk = (1.0 - bb*bb) * 14.0 * hmi_hk;
	grcrit = aa    + 0.7*(bb + 1.0);
	grc_hk = aa_hk + 0.7* bb_hk;
	gr = log10(rt);
	gr_rt = 1.0 / (2.3025851*rt);
	if(gr < grcrit-dgr) 
		
		%//----- no amplification for rtheta < rcrit
		ax    = 0.0;
		ax_hk = 0.0;
		ax_th = 0.0;
		ax_rt = 0.0;
	
    else
		
		%//----- set steep cubic ramp used to turn on ax smoothly as rtheta 
		%//-     exceeds rcrit (previously, this was done discontinuously).
		%//-     the ramp goes between  -dgr < log10(rtheta/rcrit) < dgr
		
		rnorm = (gr - (grcrit-dgr)) / (2.0*dgr);
		rn_hk =     -  grc_hk       / (2.0*dgr);
		rn_rt =  gr_rt              / (2.0*dgr);
		
		if(rnorm >= 1.0) 
			rfac    = 1.0;
			rfac_hk = 0.0;
			rfac_rt = 0.0;
		
        else
			rfac    = 3.0*rnorm*rnorm - 2.0*rnorm*rnorm*rnorm;
			rfac_rn = 6.0*rnorm    - 6.0*rnorm*rnorm;
			
			rfac_hk = rfac_rn*rn_hk;
			rfac_rt = rfac_rn*rn_rt;
        end
		
		%//----- amplification envelope slope correlation for falkner-skan
		f_arg  = 3.87*hmi    - 2.52;
		arg_hk = 3.87*hmi_hk;
		
		ex    = exp(-f_arg*f_arg);
		ex_hk = ex * (-2.0*f_arg*arg_hk);
		
		dadr    = 0.028*(hk-1.0) - 0.0345*ex;
		dadr_hk = 0.028           - 0.0345*ex_hk;
		
		%//----- new m(h) correlation    1 march 91
		af = -0.05 + 2.7*hmi -  5.5*hmi*hmi + 3.0*hmi*hmi*hmi;
		af_hmi =      2.7     - 11.0*hmi    + 9.0*hmi*hmi;
		af_hk = af_hmi*hmi_hk;
		
		ax    = (af   *dadr/th                ) * rfac;
		ax_hk = (af_hk*dadr/th + af*dadr_hk/th) * rfac...
			  + (af   *dadr/th                ) * rfac_hk;
		ax_th = -(ax)/th;
		ax_rt =  (af   *dadr/th               ) * rfac_rt;
	
    end
	
% 	return true;
% }
end
function [ ax,ax_hk1, ax_t1, ax_rt1, ax_a1,ax_hk2,ax_t2,ax_rt2,ax_a2]=axset(hk1, t1, rt1, a1,hk2, t2, rt2, a2, acrit)
% ////--------------------------------------
% //     axset function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% % {
% //----------------------------------------------------------
% //     returns average amplification ax over interval 1..2
% //----------------------------------------------------------
% //
% //==========================
% //---- 2nd-order
% 	double ax1, ax2, ax1_hk1, ax1_t1, ax1_rt1;
% 	double ax2_hk2, ax2_t2, ax2_rt2, axsq;
% 	double axa, axa_ax1, axa_ax2;
% 	double exn, exn_a1, exn_a2, dax, dax_a1, dax_a2, dax_t1, dax_t2; 
% 	double f_arg;//ex arg
	[ax1, ax1_hk1, ax1_t1, ax1_rt1]=dampl(hk1, t1, rt1);
	[ax2, ax2_hk2, ax2_t2, ax2_rt2]=dampl(hk2, t2, rt2);
% 	//---- rms-average version (seems a little better on coarse grids)
	axsq = 0.5*(ax1*ax1 + ax2*ax2);
	if(axsq <= 0.0) 
		axa = 0.0;
		axa_ax1 = 0.0;
		axa_ax2 = 0.0;
	
    else
		axa = sqrt(axsq);
		axa_ax1 = 0.5*ax1/axa;
		axa_ax2 = 0.5*ax2/axa;
    end
	
% 	//----- small additional term to ensure  dn/dx > 0  near  n = ncrit
	f_arg = min(20.0*(acrit-0.5*(a1+a2)) , 20.0);
	if(f_arg<=0.0) 
		exn    = 1.0;
		exn_a1 = 0.0;
		exn_a2 = 0.0;
	
    else
		exn    = exp(-f_arg);
		exn_a1 =  20.0*0.5*exn;
		exn_a2 =  20.0*0.5*exn;
    end
	
	dax    = exn    * 0.002/(t1+t2);
	dax_a1 = exn_a1 * 0.002/(t1+t2);
	dax_a2 = exn_a2 * 0.002/(t1+t2);
	dax_t1 = -dax/(t1+t2);
	dax_t2 = -dax/(t1+t2);
	
% 	//==========================
	
	ax     = axa             + dax;
	ax_hk1 = axa_ax1*ax1_hk1;
	ax_t1  = axa_ax1*ax1_t1  + dax_t1;
	ax_rt1 = axa_ax1*ax1_rt1;
	ax_a1  =                   dax_a1;
	
	ax_hk2 = axa_ax2*ax2_hk2;
	ax_t2  = axa_ax2*ax2_t2  + dax_t2;
	ax_rt2 = axa_ax2*ax2_rt2;
	ax_a2  =                   dax_a2;
	
% 	return true;
% }
end
function [p,blsav]=restoreblData(icom,blsav,p)
% ////--------------------------------------
% //     restoreblData function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

if (icom==1)
    p.x1     = blsav(icom).xz;
    p.u1     = blsav(icom).uz;
    p.t1     = blsav(icom).tz;
    p.d1     = blsav(icom).dz;
    p.s1     = blsav(icom).sz;
    p.ampl1  = blsav(icom).amplz;
    p.u1_uei = blsav(icom).uz_uei;
    p.u1_ms  = blsav(icom).uz_ms;
    p.dw1    = blsav(icom).dwz;
    p.h1     = blsav(icom).hz;
    p.h1_t1  = blsav(icom).hz_tz;
    p.h1_d1  = blsav(icom).hz_dz;
    p.m1     = blsav(icom).mz;
    p.m1_u1  = blsav(icom).mz_uz;
    p.m1_ms  = blsav(icom).mz_ms;
    p.r1     = blsav(icom).rz;
    p.r1_u1  = blsav(icom).rz_uz;
    p.r1_ms  = blsav(icom).rz_ms;
    p.v1     = blsav(icom).vz;
    p.v1_u1  = blsav(icom).vz_uz;
    p.v1_ms  = blsav(icom).vz_ms;
    p.v1_re  = blsav(icom).vz_re;
    p.hk1    = blsav(icom).hkz;
    p.hk1_u1 = blsav(icom).hkz_uz;
    p.hk1_t1 = blsav(icom).hkz_tz;
    p.hk1_d1 = blsav(icom).hkz_dz;
    p.hk1_ms = blsav(icom).hkz_ms;
    p.hs1    = blsav(icom).hsz;
    p.hs1_u1 = blsav(icom).hsz_uz;
    p.hs1_t1 = blsav(icom).hsz_tz;
    p.hs1_d1 = blsav(icom).hsz_dz;
    p.hs1_ms = blsav(icom).hsz_ms;
    p.hs1_re = blsav(icom).hsz_re;
    p.hc1    = blsav(icom).hcz;
    p.hc1_u1 = blsav(icom).hcz_uz;
    p.hc1_t1 = blsav(icom).hcz_tz;
    p.hc1_d1 = blsav(icom).hcz_dz;
    p.hc1_ms = blsav(icom).hcz_ms;
    p.rt1    = blsav(icom).rtz;
    p.rt1_u1 = blsav(icom).rtz_uz;
    p.rt1_t1 = blsav(icom).rtz_tz;
    p.rt1_ms = blsav(icom).rtz_ms;
    p.rt1_re = blsav(icom).rtz_re;
    p.cf1    = blsav(icom).cfz;
    p.cf1_u1 = blsav(icom).cfz_uz;
    p.cf1_t1 = blsav(icom).cfz_tz;
    p.cf1_d1 = blsav(icom).cfz_dz;
    p.cf1_ms = blsav(icom).cfz_ms;
    p.cf1_re = blsav(icom).cfz_re;
    p.di1    = blsav(icom).diz;
    p.di1_u1 = blsav(icom).diz_uz;
    p.di1_t1 = blsav(icom).diz_tz;
    p.di1_d1 = blsav(icom).diz_dz;
    p.di1_s1 = blsav(icom).diz_sz;
    p.di1_ms = blsav(icom).diz_ms;
    p.di1_re = blsav(icom).diz_re;
    p.us1    = blsav(icom).usz;
    p.us1_u1 = blsav(icom).usz_uz;
    p.us1_t1 = blsav(icom).usz_tz;
    p.us1_d1 = blsav(icom).usz_dz;
    p.us1_ms = blsav(icom).usz_ms;
    p.us1_re = blsav(icom).usz_re;
    p.cq1    = blsav(icom).cqz;
    p.cq1_u1 = blsav(icom).cqz_uz;
    p.cq1_t1 = blsav(icom).cqz_tz;
    p.cq1_d1 = blsav(icom).cqz_dz;
    p.cq1_ms = blsav(icom).cqz_ms;
    p.cq1_re = blsav(icom).cqz_re;
    p.de1    = blsav(icom).dez;
    p.de1_u1 = blsav(icom).dez_uz;
    p.de1_t1 = blsav(icom).dez_tz;
    p.de1_d1 = blsav(icom).dez_dz;
    p.de1_ms = blsav(icom).dez_ms;
end
if (icom==2)
    p.x2     = blsav(icom).xz;
    p.u2     = blsav(icom).uz;
    p.t2     = blsav(icom).tz;
    p.d2     = blsav(icom).dz;
    p.s2     = blsav(icom).sz;
    p.ampl2  = blsav(icom).amplz;
    p.u2_uei = blsav(icom).uz_uei;
    p.u2_ms  = blsav(icom).uz_ms;
    p.dw2    = blsav(icom).dwz;
    p.h2     = blsav(icom).hz;
    p.h2_t2  = blsav(icom).hz_tz;
    p.h2_d2  = blsav(icom).hz_dz;
    p.m2     = blsav(icom).mz;
    p.m2_u2  = blsav(icom).mz_uz;
    p.m2_ms  = blsav(icom).mz_ms;
    p.r2     = blsav(icom).rz;
    p.r2_u2  = blsav(icom).rz_uz;
    p.r2_ms  = blsav(icom).rz_ms;
    p.v2     = blsav(icom).vz;
    p.v2_u2  = blsav(icom).vz_uz;
    p.v2_ms  = blsav(icom).vz_ms;
    p.v2_re  = blsav(icom).vz_re;
    p.hk2    = blsav(icom).hkz;
    p.hk2_u2 = blsav(icom).hkz_uz;
    p.hk2_t2 = blsav(icom).hkz_tz;
    p.hk2_d2 = blsav(icom).hkz_dz;
    p.hk2_ms = blsav(icom).hkz_ms;
    p.hs2    = blsav(icom).hsz;
    p.hs2_u2 = blsav(icom).hsz_uz;
    p.hs2_t2 = blsav(icom).hsz_tz;
    p.hs2_d2 = blsav(icom).hsz_dz;
    p.hs2_ms = blsav(icom).hsz_ms;
    p.hs2_re = blsav(icom).hsz_re;
    p.hc2    = blsav(icom).hcz;
    p.hc2_u2 = blsav(icom).hcz_uz;
    p.hc2_t2 = blsav(icom).hcz_tz;
    p.hc2_d2 = blsav(icom).hcz_dz;
    p.hc2_ms = blsav(icom).hcz_ms;
    p.rt2    = blsav(icom).rtz;
    p.rt2_u2 = blsav(icom).rtz_uz;
    p.rt2_t2 = blsav(icom).rtz_tz;
    p.rt2_ms = blsav(icom).rtz_ms;
    p.rt2_re = blsav(icom).rtz_re;
    p.cf2    = blsav(icom).cfz;
    p.cf2_u2 = blsav(icom).cfz_uz;
    p.cf2_t2 = blsav(icom).cfz_tz;
    p.cf2_d2 = blsav(icom).cfz_dz;
    p.cf2_ms = blsav(icom).cfz_ms;
    p.cf2_re = blsav(icom).cfz_re;
    p.di2    = blsav(icom).diz;
    p.di2_u2 = blsav(icom).diz_uz;
    p.di2_t2 = blsav(icom).diz_tz;
    p.di2_d2 = blsav(icom).diz_dz;
    p.di2_s2 = blsav(icom).diz_sz;
    p.di2_ms = blsav(icom).diz_ms;
    p.di2_re = blsav(icom).diz_re;
    p.us2    = blsav(icom).usz;
    p.us2_u2 = blsav(icom).usz_uz;
    p.us2_t2 = blsav(icom).usz_tz;
    p.us2_d2 = blsav(icom).usz_dz;
    p.us2_ms = blsav(icom).usz_ms;
    p.us2_re = blsav(icom).usz_re;
    p.cq2    = blsav(icom).cqz;
    p.cq2_u2 = blsav(icom).cqz_uz;
    p.cq2_t2 = blsav(icom).cqz_tz;
    p.cq2_d2 = blsav(icom).cqz_dz;
    p.cq2_ms = blsav(icom).cqz_ms;
    p.cq2_re = blsav(icom).cqz_re;
    p.de2    = blsav(icom).dez;
    p.de2_u2 = blsav(icom).dez_uz;
    p.de2_t2 = blsav(icom).dez_tz;
    p.de2_d2 = blsav(icom).dez_dz;
    p.de2_ms = blsav(icom).dez_ms;
end
% 	return true;
end
function [blsav,p]=saveblData(icom,p,blsav)
% ////--------------------------------------
% //     saveblData function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

	if(icom==1) 
		blsav(icom).xz     = p.x1;
		blsav(icom).uz     = p.u1;
		blsav(icom).tz     = p.t1;
		blsav(icom).dz     = p.d1;
		blsav(icom).sz     = p.s1;
		blsav(icom).amplz  = p.ampl1; 
		blsav(icom).uz_uei = p.u1_uei; 
		blsav(icom).uz_ms  = p.u1_ms; 
		blsav(icom).dwz    = p.dw1; 
		blsav(icom).hz     = p.h1; 
		blsav(icom).hz_tz  = p.h1_t1; 
		blsav(icom).hz_dz  = p.h1_d1; 
		blsav(icom).mz     = p.m1; 
		blsav(icom).mz_uz  = p.m1_u1;                            
		blsav(icom).mz_ms  = p.m1_ms; 
		blsav(icom).rz     = p.r1; 
		blsav(icom).rz_uz  = p.r1_u1;                           
		blsav(icom).rz_ms  = p.r1_ms; 
		blsav(icom).vz     = p.v1; 
		blsav(icom).vz_uz  = p.v1_u1;                            
		blsav(icom).vz_ms  = p.v1_ms;  
		blsav(icom).vz_re  = p.v1_re; 
		blsav(icom).hkz    = p.hk1; 
		blsav(icom).hkz_uz = p.hk1_u1; 
		blsav(icom).hkz_tz = p.hk1_t1; 
		blsav(icom).hkz_dz = p.hk1_d1;
		blsav(icom).hkz_ms = p.hk1_ms; 
		blsav(icom).hsz    = p.hs1; 
		blsav(icom).hsz_uz = p.hs1_u1; 
		blsav(icom).hsz_tz = p.hs1_t1; 
		blsav(icom).hsz_dz = p.hs1_d1;         
		blsav(icom).hsz_ms = p.hs1_ms; 
		blsav(icom).hsz_re = p.hs1_re;
		blsav(icom).hcz    = p.hc1;
		blsav(icom).hcz_uz = p.hc1_u1;
		blsav(icom).hcz_tz = p.hc1_t1;
		blsav(icom).hcz_dz = p.hc1_d1;
		blsav(icom).hcz_ms = p.hc1_ms;
		blsav(icom).rtz    = p.rt1;
		blsav(icom).rtz_uz = p.rt1_u1;
		blsav(icom).rtz_tz = p.rt1_t1;
		blsav(icom).rtz_ms = p.rt1_ms;
		blsav(icom).rtz_re = p.rt1_re;
		blsav(icom).cfz    = p.cf1;
		blsav(icom).cfz_uz = p.cf1_u1;
		blsav(icom).cfz_tz = p.cf1_t1;
		blsav(icom).cfz_dz = p.cf1_d1;
		blsav(icom).cfz_ms = p.cf1_ms;
		blsav(icom).cfz_re = p.cf1_re;
		blsav(icom).diz    = p.di1;
		blsav(icom).diz_uz = p.di1_u1;
		blsav(icom).diz_tz = p.di1_t1;
		blsav(icom).diz_dz = p.di1_d1;
		blsav(icom).diz_sz = p.di1_s1;
		blsav(icom).diz_ms = p.di1_ms;
		blsav(icom).diz_re = p.di1_re;
		blsav(icom).usz    = p.us1;
		blsav(icom).usz_uz = p.us1_u1;
		blsav(icom).usz_tz = p.us1_t1;
		blsav(icom).usz_dz = p.us1_d1;
		blsav(icom).usz_ms = p.us1_ms;
		blsav(icom).usz_re = p.us1_re; 
		blsav(icom).cqz    = p.cq1; 
		blsav(icom).cqz_uz = p.cq1_u1; 
		blsav(icom).cqz_tz = p.cq1_t1; 
		blsav(icom).cqz_dz = p.cq1_d1;         
		blsav(icom).cqz_ms = p.cq1_ms;
		blsav(icom).cqz_re = p.cq1_re;
		blsav(icom).dez    = p.de1; 
		blsav(icom).dez_uz = p.de1_u1; 
		blsav(icom).dez_tz = p.de1_t1;
		blsav(icom).dez_dz = p.de1_d1;         
		blsav(icom).dez_ms = p.de1_ms;
    
    else
		blsav(icom).xz     = p.x2;
		blsav(icom).uz     = p.u2;
		blsav(icom).tz     = p.t2;
		blsav(icom).dz     = p.d2;
		blsav(icom).sz     = p.s2;
		blsav(icom).amplz  = p.ampl2; 
		blsav(icom).uz_uei = p.u2_uei; 
		blsav(icom).uz_ms  = p.u2_ms; 
		blsav(icom).dwz    = p.dw2; 
		blsav(icom).hz     = p.h2; 
		blsav(icom).hz_tz  = p.h2_t2; 
		blsav(icom).hz_dz  = p.h2_d2; 
		blsav(icom).mz     = p.m2; 
		blsav(icom).mz_uz  = p.m2_u2;                            
		blsav(icom).mz_ms  = p.m2_ms; 
		blsav(icom).rz     = p.r2; 
		blsav(icom).rz_uz  = p.r2_u2;                           
		blsav(icom).rz_ms  = p.r2_ms; 
		blsav(icom).vz     = p.v2; 
		blsav(icom).vz_uz  = p.v2_u2;                            
		blsav(icom).vz_ms  = p.v2_ms;  
		blsav(icom).vz_re  = p.v2_re; 
		blsav(icom).hkz    = p.hk2; 
		blsav(icom).hkz_uz = p.hk2_u2; 
		blsav(icom).hkz_tz = p.hk2_t2; 
		blsav(icom).hkz_dz = p.hk2_d2;
		blsav(icom).hkz_ms = p.hk2_ms; 
		blsav(icom).hsz    = p.hs2; 
		blsav(icom).hsz_uz = p.hs2_u2; 
		blsav(icom).hsz_tz = p.hs2_t2; 
		blsav(icom).hsz_dz = p.hs2_d2;         
		blsav(icom).hsz_ms = p.hs2_ms; 
		blsav(icom).hsz_re = p.hs2_re;
		blsav(icom).hcz    = p.hc2;
		blsav(icom).hcz_uz = p.hc2_u2;
		blsav(icom).hcz_tz = p.hc2_t2;
		blsav(icom).hcz_dz = p.hc2_d2;
		blsav(icom).hcz_ms = p.hc2_ms;
		blsav(icom).rtz    = p.rt2;
		blsav(icom).rtz_uz = p.rt2_u2;
		blsav(icom).rtz_tz = p.rt2_t2;
		blsav(icom).rtz_ms = p.rt2_ms;
		blsav(icom).rtz_re = p.rt2_re;
		blsav(icom).cfz    = p.cf2;
		blsav(icom).cfz_uz = p.cf2_u2;
		blsav(icom).cfz_tz = p.cf2_t2;
		blsav(icom).cfz_dz = p.cf2_d2;
		blsav(icom).cfz_ms = p.cf2_ms;
		blsav(icom).cfz_re = p.cf2_re;
		blsav(icom).diz    = p.di2;
		blsav(icom).diz_uz = p.di2_u2;
		blsav(icom).diz_tz = p.di2_t2;
		blsav(icom).diz_dz = p.di2_d2;
		blsav(icom).diz_sz = p.di2_s2;
		blsav(icom).diz_ms = p.di2_ms;
		blsav(icom).diz_re = p.di2_re;
		blsav(icom).usz    = p.us2;
		blsav(icom).usz_uz = p.us2_u2;
		blsav(icom).usz_tz = p.us2_t2;
		blsav(icom).usz_dz = p.us2_d2;
		blsav(icom).usz_ms = p.us2_ms;
		blsav(icom).usz_re = p.us2_re; 
		blsav(icom).cqz    = p.cq2; 
		blsav(icom).cqz_uz = p.cq2_u2; 
		blsav(icom).cqz_tz = p.cq2_t2; 
		blsav(icom).cqz_dz = p.cq2_d2;         
		blsav(icom).cqz_ms = p.cq2_ms;
		blsav(icom).cqz_re = p.cq2_re;
		blsav(icom).dez    = p.de2; 
		blsav(icom).dez_uz = p.de2_u2; 
		blsav(icom).dez_tz = p.de2_t2;
		blsav(icom).dez_dz = p.de2_d2;         
		blsav(icom).dez_ms = p.de2_ms;
    end
% 	return true;
% }e
end
function p=trchek(p)
% ////--------------------------------------
% //     trchek function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% //----------------------------------------------------------------
% //     new second-order version:  december 1994.
% //
% //     checks if transition occurs in the current interval p.x1..p.x2.
% //     if transition occurs, then set transition location p.xt, and 
% //     its sensitivities to "1" and "2" variables.  if no transition, 
% //     set amplification p.ampl2.
% //
% //     solves the implicit amplification equation for n2:
% //
% //       n2 - n1     n'(p.xt,nt) + n'(p.x1,n1)
% //       -------  =  ---------------------
% //       p.x2 - p.x1               2
% //
% //     in effect, a 2-point central difference is used between
% //     p.x1..p.x2 (no transition), or p.x1..p.xt (transition).  the switch
% //     is done by defining p.xt,nt in the equation above depending
% //     on whether n2 exceeds ncrit.
% //
% //  if n2<ncrit:  nt=n2    , p.xt=p.x2                  (no transition)
% //
% //  if n2>ncrit:  nt=ncrit , p.xt=(ncrit-n1)/(n2-n1)  (transition)
% //
% //----------------------------------------------------------------
% 	CString str;

	 daeps =0.00005;
% 	int itam;
% 	 p.ax_hk1, p.ax_t1, p.ax_a1, p.ax_hk2, p.ax_t2, p.ax_rt2, p.ax_a2;
% 	 p.amplt, p.sfa, p.sfa_a1, p.sfa_a2, p.sfx;
% 	 p.sfx_x1, p.sfx_x2, p.sfx_xf;
% 	 p.tt, p.dt, p.ut,  p.amsave;
% 	 p.ax, p.ax_rt1,res, res_a2;
% 	 p.da2, p.dxt, p.tt_t1, p.dt_d1, p.ut_u1;
% 	 p.tt_t2, p.dt_d2, p.ut_u2, p.tt_a1, p.dt_a1;
% 	 p.ut_a1, p.tt_x1, p.dt_x1, p.ut_x1, p.tt_x2, p.dt_x2, p.ut_x2, p.tt_xf, p.dt_xf, p.ut_xf;
% 	 p.ax_d1, p.ax_u1, p.ax_x1, p.ax_d2, p.ax_u2, p.ax_x2, p.ax_xf, p.ax_ms, p.ax_re;
% 	 p.z_ax, p.z_a1, p.z_t1, p.z_d1, p.z_u1, p.z_x1, p.z_a2, p.z_t2, p.z_d2, p.z_u2,
% 		  p.z_x2, p.z_xf, p.z_ms, p.z_re;
% 
% 	//intialization to 0.0 added arcds to avoid level 4 warnings at compile time
	 p.ax_at = 0.0;
p.ax_rtt = 0.0;
p.ax_tt = 0.0;
p.ax_hkt = 0.0;
p.amplt_a2 = 0.0; 
p.wf1 = 0.0;
p.wf1_a1 = 0.0;
p.wf1_a2 = 0.0;
p.wf1_xf = 0.0;
p.wf1_x1 = 0.0;
p.wf1_x2 = 0.0;
p.wf2 = 0.0;
p.wf2_a1 = 0.0;
p.wf2_a2 = 0.0;
p.wf2_xf = 0.0;
p.wf2_x1 = 0.0;
p.wf2_x2 = 0.0;
p.xt_a2 = 0.0;
p.dt_a2 = 0.0;
p.tt_a2 = 0.0;
p.ut_a2 = 0.0;
p.hkt =0.0 ;
p.hkt_tt = 0.0;
p.hkt_dt = 0.0;
p.hkt_ut = 0.0;
p.hkt_ms = 0.0;
p.rtt_tt = 0.0;
p.rtt_ut = 0.0;
p.rtt_ms = 0.0;
p.rtt= 0.0;
p.rtt_re = 0.0;
% 	//---- save variables and sensitivities at ibl ("2") for future restoration
    blsav(2).d=0;
	%//---- save variables and sensitivities at ibl ("2") for future restoration
 	[blsav,p]=saveblData(2,p,blsav);
	
	%//---- calculate average amplification rate p.ax over p.x1..p.x2 interval
%     	axset(hk1, p.t1, rt1, p.ampl1, p.hk2, p.t2, p.rt2, p.ampl2, p.amcrit,
% 		p.ax, p.ax_hk1, p.ax_t1, p.ax_rt1, p.ax_a1, p.ax_hk2, p.ax_t2, p.ax_rt2, p.ax_a2 );

	[p.ax, p.ax_hk1, p.ax_t1, p.ax_rt1, p.ax_a1, p.ax_hk2, p.ax_t2, p.ax_rt2, p.ax_a2]=axset(p.hk1, p.t1, p.rt1, p.ampl1, p.hk2, p.t2, p.rt2, p.ampl2, p.amcrit );
	
	%//---- set initial guess for iterate n2 (p.ampl2) at p.x2
	p.ampl2 = p.ampl1 + p.ax*(p.x2-p.x1);
	
% 	//---- solve implicit system for amplification p.ampl2
	for(itam=1:30)
		
% 		//---- define weighting factors p.wf1,p.wf2 for defining "t" quantities from 1,2
		if(p.ampl2 <= p.amcrit) 
% 			//------ there is no transition yet,  "t" is the same as "2"
			p.amplt    = p.ampl2;
			p.amplt_a2 = 1.0;
			p.sfa    = 1.0;
			p.sfa_a1 = 0.0;
			p.sfa_a2 = 0.0;
		
        else
% 			//------ there is transition in p.x1..p.x2, "t" is set from n1, n2
			p.amplt    = p.amcrit;
			p.amplt_a2 = 0.0;
			p.sfa    = (p.amplt - p.ampl1)/(p.ampl2-p.ampl1);
			p.sfa_a1 = ( p.sfa  - 1.0  )/(p.ampl2-p.ampl1);
			p.sfa_a2 = (      - p.sfa  )/(p.ampl2-p.ampl1);
        end
		
		if(p.xiforc<p.x2) 
			p.sfx    = (p.xiforc - p.x1 )/(p.x2-p.x1);
			p.sfx_x1 = (p.sfx    - 1.0)/(p.x2-p.x1);
			p.sfx_x2 = (       - p.sfx)/(p.x2-p.x1);
			p.sfx_xf =  1.0          /(p.x2-p.x1);
		
        else
			p.sfx    = 1.0;
			p.sfx_x1 = 0.0;
			p.sfx_x2 = 0.0;
			p.sfx_xf = 0.0;
        end
		
% 		//---- set weighting factor from free or forced transition
		if(p.sfa<p.sfx) 
			p.wf2    = p.sfa;
			p.wf2_a1 = p.sfa_a1;
			p.wf2_a2 = p.sfa_a2;
			p.wf2_x1 = 0.0;
			p.wf2_x2 = 0.0;
			p.wf2_xf = 0.0;
		
        else
			p.wf2    = p.sfx;
			p.wf2_a1 = 0.0;
			p.wf2_a2 = 0.0;
			p.wf2_x1 = p.sfx_x1;
			p.wf2_x2 = p.sfx_x2;
			p.wf2_xf = p.sfx_xf;
        end
		
		p.wf1    = 1.0 - p.wf2;
		p.wf1_a1 =     - p.wf2_a1;
		p.wf1_a2 =     - p.wf2_a2;
		p.wf1_x1 =     - p.wf2_x1;
		p.wf1_x2 =     - p.wf2_x2;
		p.wf1_xf =     - p.wf2_xf;
		
% 		//---- interpolate bl variables to p.xt
		p.xt    = p.x1*p.wf1    + p.x2*p.wf2;
		p.tt    = p.t1*p.wf1    + p.t2*p.wf2;
		p.dt    = p.d1*p.wf1    + p.d2*p.wf2;
		p.ut    = p.u1*p.wf1    + p.u2*p.wf2;
		
		p.xt_a2 = p.x1*p.wf1_a2 + p.x2*p.wf2_a2;
		p.tt_a2 = p.t1*p.wf1_a2 + p.t2*p.wf2_a2;
		p.dt_a2 = p.d1*p.wf1_a2 + p.d2*p.wf2_a2;
		p.ut_a2 = p.u1*p.wf1_a2 + p.u2*p.wf2_a2;
		
% 		//---- temporarily set "2" variables from "t" for blkin
		p.x2 = p.xt;
		p.t2 = p.tt;
		p.d2 = p.dt;
		p.u2 = p.ut;
		
% 		//---- calculate laminar secondary "t" variables p.hkt, p.rtt
		p=blkin(p);
	
		p.hkt    = p.hk2;
		p.hkt_tt = p.hk2_t2;
		p.hkt_dt = p.hk2_d2;
		p.hkt_ut = p.hk2_u2;
		p.hkt_ms = p.hk2_ms;
		
		p.rtt    = p.rt2;
		p.rtt_tt = p.rt2_t2;
		p.rtt_ut = p.rt2_u2;
		p.rtt_ms = p.rt2_ms;
		p.rtt_re = p.rt2_re;
		
% 		//---- restore clobbered "2" variables, except for p.ampl2
		p.amsave = p.ampl2;

 		[p,blsav]=restoreblData(2,blsav,p);
		
		p.ampl2 = p.amsave;
	
		%//---- calculate amplification rate p.ax over current x1-xt interval
		[p.ax, p.ax_hk1, p.ax_t1, p.ax_rt1, p.ax_a1,p.ax_hkt, p.ax_tt, p.ax_rtt, p.ax_at]=axset(p.hk1, p.t1, p.rt1, p.ampl1, p.hkt, p.tt, p.rtt, p.amplt,p.amcrit);


		
% 		//---- punch out early if there is no amplification here
		if(p.ax <= 0.0) 
            break
        end
		
% 		//---- set sensitivity of p.ax(a2)
		p.ax_a2 = (p.ax_hkt*p.hkt_tt + p.ax_tt + p.ax_rtt*p.rtt_tt)*p.tt_a2...
			+ (p.ax_hkt*p.hkt_dt                        )*p.dt_a2...
			+ (p.ax_hkt*p.hkt_ut         + p.ax_rtt*p.rtt_ut)*p.ut_a2...
			+  p.ax_at                                 *p.amplt_a2;
		
% 		//---- residual for implicit p.ampl2 definition (amplification equation)
		res    = p.ampl2 - p.ampl1 - p.ax   *(p.x2-p.x1);
		res_a2 = 1.0          - p.ax_a2*(p.x2-p.x1);
		
		p.da2 = -res/res_a2;
		
		p.rlx = 1.0;
		p.dxt = p.xt_a2*p.da2;
		
		if(p.rlx*abs(p.dxt/(p.x2-p.x1)) > 0.05) 
            p.rlx = 0.05*abs((p.x2-p.x1)/p.dxt);
        end
		if(p.rlx*abs(p.da2)         > 1.0)  
            p.rlx = 1.0 *abs(   1.0 /p.da2);
        end
		
% 		//---- check if converged
		if(abs(p.da2) < daeps) 
            break
        end

		if((p.ampl2>p.amcrit && p.ampl2+p.rlx*p.da2<p.amcrit)||	(p.ampl2<p.amcrit && p.ampl2+p.rlx*p.da2>p.amcrit)    ) 
% 			//------ limited newton step so p.ampl2 doesn't step across p.amcrit either way
			p.ampl2 = p.amcrit;
		else
% 			//------ regular newton step
			p.ampl2 = p.ampl2 + p.rlx*p.da2;
        end
    end
   if itam==30
	disp('trchek2 - n2 convergence failed');
% 	str.Format("trchek2 - n2 convergence failed\r\n");
% 	if(m_bTrace)pXFile->WriteString(str);
   end

% stop101:

% 	//---- test for free or forced transition
	p.trfree = (p.ampl2 >= p.amcrit);
	p.trforc = (p.xiforc>p.x1) && (p.xiforc<=p.x2);

% 	//---- set transition interval flag
	p.tran = (p.trforc || p.trfree);

	if(p.tran==false) 
        return;
    end

% 	//---- resolve if both forced and free transition
	if(p.trfree && p.trforc) 
		p.trforc = p.xiforc < p.xt;
		p.trfree = p.xiforc >= p.xt;
    end

	if(p.trforc) 
% 		//----- if forced transition, then p.xt is prescribed,
% 		//-     no sense calculating the sensitivities, since we know them...
		p.xt = p.xiforc;
		p.xt_a1 = 0.0;
		p.xt_x1 = 0.0;
		p.xt_t1 = 0.0;
		p.xt_d1 = 0.0;
		p.xt_u1 = 0.0;
		p.xt_x2 = 0.0;
		p.xt_t2 = 0.0;
		p.xt_d2 = 0.0;
		p.xt_u2 = 0.0;
		p.xt_ms = 0.0;
		p.xt_re = 0.0;
		p.xt_xf = 1.0;
% 		return true;
    end

% 	//---- free transition ... set sensitivities of p.xt

	p.xt_x1 =    p.wf1;
	p.tt_t1 =    p.wf1;
	p.dt_d1 =    p.wf1;
	p.ut_u1 =    p.wf1;

	p.xt_x2 =                p.wf2;
	p.tt_t2 =                p.wf2;
	p.dt_d2 =                p.wf2;
	p.ut_u2 =                p.wf2;

	p.xt_a1 = p.x1*p.wf1_a1 + p.x2*p.wf2_a1;
	p.tt_a1 = p.t1*p.wf1_a1 + p.t2*p.wf2_a1;
	p.dt_a1 = p.d1*p.wf1_a1 + p.d2*p.wf2_a1;
	p.ut_a1 = p.u1*p.wf1_a1 + p.u2*p.wf2_a1;

	p.xt_x1 = p.x1*p.wf1_x1 + p.x2*p.wf2_x1 + p.xt_x1;
	p.tt_x1 = p.t1*p.wf1_x1 + p.t2*p.wf2_x1;
	p.dt_x1 = p.d1*p.wf1_x1 + p.d2*p.wf2_x1;
	p.ut_x1 = p.u1*p.wf1_x1 + p.u2*p.wf2_x1;

	p.xt_x2 = p.x1*p.wf1_x2 + p.x2*p.wf2_x2 + p.xt_x2;
	p.tt_x2 = p.t1*p.wf1_x2 + p.t2*p.wf2_x2;
	p.dt_x2 = p.d1*p.wf1_x2 + p.d2*p.wf2_x2;
	p.ut_x2 = p.u1*p.wf1_x2 + p.u2*p.wf2_x2;

	p.xt_xf = p.x1*p.wf1_xf + p.x2*p.wf2_xf;
	p.tt_xf = p.t1*p.wf1_xf + p.t2*p.wf2_xf;
	p.dt_xf = p.d1*p.wf1_xf + p.d2*p.wf2_xf;
	p.ut_xf = p.u1*p.wf1_xf + p.u2*p.wf2_xf;

% 	//---- at this point, p.ax = p.ax( hk1, p.t1, rt1, a1, p.hkt, p.tt, p.rtt, at )

% 	//---- set sensitivities of p.ax( p.t1 p.d1 p.u1 a1 p.t2 p.d2 p.u2 a2 ms re )
	p.ax_t1 =  p.ax_hk1*p.hk1_t1 + p.ax_t1 + p.ax_rt1*p.rt1_t1...
	+ (p.ax_hkt*p.hkt_tt + p.ax_tt + p.ax_rtt*p.rtt_tt)*p.tt_t1;
	p.ax_d1 =  p.ax_hk1*p.hk1_d1...
	+ (p.ax_hkt*p.hkt_dt                        )*p.dt_d1;
	p.ax_u1 =  p.ax_hk1*p.hk1_u1         + p.ax_rt1*p.rt1_u1...
	+ (p.ax_hkt*p.hkt_ut         + p.ax_rtt*p.rtt_ut)*p.ut_u1;
	p.ax_a1 =  p.ax_a1...
	+ (p.ax_hkt*p.hkt_tt + p.ax_tt + p.ax_rtt*p.rtt_tt)*p.tt_a1...
	+ (p.ax_hkt*p.hkt_dt                        )*p.dt_a1...
	+ (p.ax_hkt*p.hkt_ut         + p.ax_rtt*p.rtt_ut)*p.ut_a1;
	p.ax_x1 = (p.ax_hkt*p.hkt_tt + p.ax_tt + p.ax_rtt*p.rtt_tt)*p.tt_x1...
	+ (p.ax_hkt*p.hkt_dt                        )*p.dt_x1...
	+ (p.ax_hkt*p.hkt_ut         + p.ax_rtt*p.rtt_ut)*p.ut_x1;

	p.ax_t2 = (p.ax_hkt*p.hkt_tt + p.ax_tt + p.ax_rtt*p.rtt_tt)*p.tt_t2;
	p.ax_d2 = (p.ax_hkt*p.hkt_dt                        )*p.dt_d2;
	p.ax_u2 = (p.ax_hkt*p.hkt_ut         + p.ax_rtt*p.rtt_ut)*p.ut_u2;
	p.ax_a2 =  p.ax_at                                 *p.amplt_a2...
	+ (p.ax_hkt*p.hkt_tt + p.ax_tt + p.ax_rtt*p.rtt_tt)*p.tt_a2...
	+ (p.ax_hkt*p.hkt_dt                        )*p.dt_a2...
	+ (p.ax_hkt*p.hkt_ut         + p.ax_rtt*p.rtt_ut)*p.ut_a2;
	p.ax_x2 = (p.ax_hkt*p.hkt_tt + p.ax_tt + p.ax_rtt*p.rtt_tt)*p.tt_x2...
	+ (p.ax_hkt*p.hkt_dt                        )*p.dt_x2...
	+ (p.ax_hkt*p.hkt_ut         + p.ax_rtt*p.rtt_ut)*p.ut_x2;

	p.ax_xf = (p.ax_hkt*p.hkt_tt + p.ax_tt + p.ax_rtt*p.rtt_tt)*p.tt_xf...
	+ (p.ax_hkt*p.hkt_dt                        )*p.dt_xf...
	+ (p.ax_hkt*p.hkt_ut         + p.ax_rtt*p.rtt_ut)*p.ut_xf;

	p.ax_ms =  p.ax_hkt*p.hkt_ms         + p.ax_rtt*p.rtt_ms...
	+  p.ax_hk1*p.hk1_ms         + p.ax_rt1*p.rt1_ms;
	p.ax_re =                          p.ax_rtt*p.rtt_re...
	+ p.ax_rt1*p.rt1_re;

% 	//---- set sensitivities of residual res
% 	//c   res  = p.ampl2 - p.ampl1 - p.ax*(p.x2-p.x1)
	p.z_ax =               -    (p.x2-p.x1);

	p.z_a1 = p.z_ax*p.ax_a1 - 1.0;
	p.z_t1 = p.z_ax*p.ax_t1;
	p.z_d1 = p.z_ax*p.ax_d1;
	p.z_u1 = p.z_ax*p.ax_u1;
	p.z_x1 = p.z_ax*p.ax_x1 + p.ax;

	p.z_a2 = p.z_ax*p.ax_a2 + 1.0;
	p.z_t2 = p.z_ax*p.ax_t2;
	p.z_d2 = p.z_ax*p.ax_d2;
	p.z_u2 = p.z_ax*p.ax_u2;
	p.z_x2 = p.z_ax*p.ax_x2 - p.ax;

	p.z_xf = p.z_ax*p.ax_xf;
	p.z_ms = p.z_ax*p.ax_ms;
	p.z_re = p.z_ax*p.ax_re;

% 	//---- set sensitivities of p.xt, with res being stationary for a2 constraint
	p.xt_a1 = p.xt_a1 - (p.xt_a2/p.z_a2)*p.z_a1;
	p.xt_t1 =       - (p.xt_a2/p.z_a2)*p.z_t1;
	p.xt_d1 =       - (p.xt_a2/p.z_a2)*p.z_d1;
	p.xt_u1 =       - (p.xt_a2/p.z_a2)*p.z_u1;
	p.xt_x1 = p.xt_x1 - (p.xt_a2/p.z_a2)*p.z_x1;
	p.xt_t2 =       - (p.xt_a2/p.z_a2)*p.z_t2;
	p.xt_d2 =       - (p.xt_a2/p.z_a2)*p.z_d2;
	p.xt_u2 =       - (p.xt_a2/p.z_a2)*p.z_u2;
	p.xt_x2 = p.xt_x2 - (p.xt_a2/p.z_a2)*p.z_x2;
	p.xt_ms =       - (p.xt_a2/p.z_a2)*p.z_ms;
	p.xt_re =       - (p.xt_a2/p.z_a2)*p.z_re;
	p.xt_xf = 0.0;

% 	return true ;
% }
end

function dstr=dslim( dstr,  thet,  msq,  hklim)
% ////--------------------------------------
% //     dslim function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------


	 h = (dstr)/thet;
% 	 hk, hk_h, hk_m;
	
	[hk, hk_h, hk_m]=hkin(h, msq);
	
	 dh = max(0.0 , hklim-hk ) / hk_h;
	dstr = (dstr) + dh*thet;
	
% 	return true;
% }
end
function r=Gauss(nn, z, r)
% ////--------------------------------------
% //     Gauss function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% /*******************************************************
%  *                                                     *
%  *   solves general nxn system in nn unknowns          *
%  *    with arbitrary number (nrhs) of righthand sides. *
%  *   assumes system is invertible...                   *
%  *    ...if it isn't, a divide by zero will result.    *
%  *                                                     *
%  *   z is the coefficient matrix...                    *
%  *     ...destroyed during solution process.           *
%  *   r is the righthand side(s)...                     *
%  *     ...replaced by the solution vector(s).          *
%  *                                                     *
%  *                              mark drela  1984       *
%  *******************************************************/
% %// arcds : only one rhs is enough ! nrhs = 1
% %// dimension z(nsiz,nsiz), r(nsiz,nrhs)
% 
% 	int loc;
% 	int np, nnpp, nt, nx, k;
% 
% 	double temp, ztmp, pivot;
% 	
	for (np=1:nn-1)
		nnpp = np+1;
		%//------ find max pivot index nx
		nx = np;
		for (nt =nnpp:nn)
			if (abs(z(nt,np))>abs(z(nx,np))) 
                nx = nt;
            end
        end

		pivot = 1.0./z(nx,np);
		
		%//------ switch pivots
		z(nx,np) = z(np,np);
		
		%//------ switch rows & normalize pivot row
		for (loc = nnpp:nn)
			temp = z(nx,loc)*pivot;
			z(nx,loc) = z(np,loc);
			z(np,loc) = temp;
        end

		temp = r(nx)*pivot;
		r(nx) = r(np);
		r(np) = temp;

		%//------ forward eliminate everything
		for (k = nnpp:nn)
			ztmp = z(k,np);
			for (loc=nnpp:nn) 
                z(k,loc) = z(k,loc) - ztmp*z(np,loc);
            end
			r(k) = r(k) - ztmp*r(np);
        end
    end
	
	%//---- solve for last row
	r(nn) = r(nn)/z(nn,nn);
	
	%//---- back substitute everything
	for (np=nn-1:-1: 1)
		nnpp = np+1;
		for(k=nnpp:nn)
			r(np) = r(np) - z(np,k)*r(k);
        end
    end
% 
% 	return true;
% }
end
function [p]=bldif(ityp,p)
% ////--------------------------------------
% //     bldif function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% {
	%//-----------------------------------------------------------
	%//     sets up the newton system coefficients and residuals
	%//
	%//        ityp = 0 :  similarity station
	%//        ityp = 1 :  laminar interval
	%//        ityp = 2 :  turp.bulent interval
	%//        ityp = 3 :  wake interval
	%//
	%//     this routine knows nothing about a transition interval,
	%//     which is taken care of by trdif.
	%//-----------------------------------------------------------
% 	int k,l;
% 	double p.hupwt,p.hdcon, p.hl, p.hd_hk1,hd_hk2, p.hlsq, p.ehh;
% 	double p.upw, p.upw_hl, p.upw_hd, p.upw_hk1, p.upw_hk2, upw_u1, p.upw_t1, p.upw_d1;
% 	double p.upw_u2, p.upw_t2, p.upw_d2, p.upw_ms;
% 	double p.dxi, p.slog, p.scc, scc_usa, p.scc_us1, scc_us2;
% 	double ax, ax_t1, ax_rp.t1, p.ax_a1, ax_hk2, ax_t2, ax_rp.t2, ax_a2;
% 	double p.rezc, p.z_ax, p.hr, p.hr_hka, p.hr_rta;
% 	double p.hl_hk1, hl_hk2, p.ax_hk1, p.sa, p.cqa, p.cfa, p.hka, p.usa, rta, p.dea, p.da, p.ald;
% 	double gcc, p.hkc, p.hkc_hka, p.hkc_rta, rezt, rezh;
% 	double btmp, hwa, ha, ma, xa, ta, p.xlog, p.ulog, p.tlog, p.hlog, p.ddlog;
% 	double z_cfx, z_ha, z_hwa, z_ma, z_xl, z_tl, z_cfm, p.z_t1, p.z_t2;
% 	double z_dix, z_hca, z_hl, z_hs1, z_hs2, z_di1, z_di2;
% 	double p.z_cfa, p.z_hka, p.z_da, p.z_sl, p.z_ul, p.z_dxi, p.z_usa, p.z_cqa, p.z_sa, p.z_dea;
% 	double p.z_upw, p.z_de1, p.z_de2, p.z_us1, p.z_us2, p.z_d1, p.z_d2, z_u1, z_u2, p.z_x1, z_x2;
% 	double p.z_s1, p.z_s2, p.z_cq1, p.z_cq2, p.z_cf1, p.z_cf2, p.z_hk1, z_hk2;
% 	double cfx, cfx_xa, cfx_ta, cfx_x1, cfx_x2, cfx_t1, cfx_t2, cfx_cf1, cfx_cf2, cfx_cfm;
% 	double xop.t1, xop.t2, hca, hsa, dix, dix_upw, cfx_upw;
% 	double p.uq, p.uq_hka, p.uq_rta, p.uq_cfa, p.uq_da, p.uq_upw, uq_t1, uq_t2;
% 	double uq_d1, uq_d2, uq_u1, uq_u2, uq_ms, uq_re;
% 	double p.f_arg;%// ex arg

	if(ityp==0)
		%//----- similarity logarithmic differences  (prescribed)
		p.xlog = 1.0;
		p.ulog = p.bule;
		p.tlog = 0.5*(1.0 - p.bule);
		p.hlog = 0.0;
		p.ddlog = 0.0;
	
    else
		%//----- usual logarithmic differences
		p.xlog = log(p.x2/p.x1);
		p.ulog = log(p.u2/p.u1);
		p.tlog = log(p.t2/p.t1);
		p.hlog = log(p.hs2/p.hs1);
		p.ddlog = 1.0;
    end
	
	for (k=1:4)
		p.vsrez(k) = 0.0;
		p.vsm(k) = 0.0;
		p.vsr(k) = 0.0;
		p.vsx(k) = 0.0;
		for (l=1:5)
			p.vs1(k,l) = 0.0; 
			p.vs2(k,l) = 0.0;
        end
    end
	
	
	%//---- set triggering constant for local upwinding
	p.hupwt = 1.0;
	
	p.hdcon  =  5.0*p.hupwt/p.hk2/p.hk2;
	p.hd_hk1 =  0.0;
	p.hd_hk2 = -p.hdcon*2.0/p.hk2;
	
	%//---- use less upwinding in the wake
	if(ityp==3) 
		p.hdcon  =  p.hupwt/p.hk2/p.hk2;
		p.hd_hk1 =  0.0;
		p.hd_hk2 = -p.hdcon*2.0/p.hk2;
	end
	%//
	%//---- local upwinding is based on local change in  log(hk-1)
	%//-    (mainly kicks in at transition)
	p.f_arg = abs((p.hk2-1.0)/(p.hk1-1.0));
	p.hl = log(p.f_arg);
	p.hl_hk1 = -1.0/(p.hk1-1.0);
	p.hl_hk2 =  1.0/(p.hk2-1.0);
	
	p.hlsq = min(p.hl*p.hl, 15.0);
	p.ehh = exp(-p.hlsq*p.hdcon);
	p.upw = 1.0 - 0.5*p.ehh;
	p.upw_hl =        p.ehh * p.hl  *p.hdcon;
	p.upw_hd =    0.5*p.ehh * p.hlsq;
	
	p.upw_hk1 = p.upw_hl*p.hl_hk1 + p.upw_hd*p.hd_hk1;
	p.upw_hk2 = p.upw_hl*p.hl_hk2 + p.upw_hd*p.hd_hk2;
	
	p.upw_u1 = p.upw_hk1*p.hk1_u1;
	p.upw_t1 = p.upw_hk1*p.hk1_t1;
	p.upw_d1 = p.upw_hk1*p.hk1_d1;
	p.upw_u2 = p.upw_hk2*p.hk2_u2;
	p.upw_t2 = p.upw_hk2*p.hk2_t2;
	p.upw_d2 = p.upw_hk2*p.hk2_d2;
	p.upw_ms = p.upw_hk1*p.hk1_ms + p.upw_hk2*p.hk2_ms;
	
	
	if(ityp==0) 
		
		%//***** le point -->  set zero amplification factor
		p.vs2(1,1) = 1.0;
		p.vsr(1)   = 0.0;
		p.vsrez(1) = -p.ampl2;
	
	else 
		if(ityp==1) 
		%//***** laminar part -->  set amplification equation
		%//----- set average amplification ax over interval p.x1..p.x2
		
		[p.ax,p.ax_hk1, p.ax_t1, p.ax_rt1, p.ax_a1,  p.ax_hk2, p.ax_t2, p.ax_rt2, p.ax_a2]=axset(p.hk1, p.t1, p.rt1, p.ampl1, ...
			  p.hk2, p.t2, p.rt2, p.ampl2,...
			  p.amcrit);
		
		p.rezc = p.ampl2 - p.ampl1 - p.ax*(p.x2-p.x1);
		p.z_ax = -(p.x2-p.x1);
		
		p.vs1(1,1) = p.z_ax* p.ax_a1  -  1.0;
		p.vs1(1,2) = p.z_ax*(p.ax_hk1*p.hk1_t1 + p.ax_t1 + p.ax_rt1*p.rt1_t1);
		p.vs1(1,3) = p.z_ax*(p.ax_hk1*p.hk1_d1                        );
		p.vs1(1,4) = p.z_ax*(p.ax_hk1*p.hk1_u1         + p.ax_rt1*p.rt1_u1);
		p.vs1(1,5) =  p.ax;
		p.vs2(1,1) = p.z_ax* p.ax_a2  +  1.0;
		p.vs2(1,2) = p.z_ax*(p.ax_hk2*p.hk2_t2 + p.ax_t2 + p.ax_rt2*p.rt2_t2);
		p.vs2(1,3) = p.z_ax*(p.ax_hk2*p.hk2_d2                        ) ;        
		p.vs2(1,4) = p.z_ax*(p.ax_hk2*p.hk2_u2         + p.ax_rt2*p.rt2_u2);
		p.vs2(1,5) = -p.ax;
		p.vsm(1)   = p.z_ax*(p.ax_hk1*p.hk1_ms         + p.ax_rt1*p.rt1_ms...
			+ p.ax_hk2*p.hk2_ms         + p.ax_rt2*p.rt2_ms);
		p.vsr(1)   = p.z_ax*(                        p.ax_rt1*p.rt1_re...
			+ p.ax_rt2*p.rt2_re);
		p.vsx(1)   = 0.0;
		p.vsrez(1) = -p.rezc;
	
        else
		
		%//***** turp.bulent part -->  set shear lag equation
		
		p.sa  = (1.0-p.upw)*p.s1  + p.upw*p.s2;
		p.cqa = (1.0-p.upw)*p.cq1 + p.upw*p.cq2;
		p.cfa = (1.0-p.upw)*p.cf1 + p.upw*p.cf2;
		p.hka = (1.0-p.upw)*p.hk1 + p.upw*p.hk2;
		
		p.usa = 0.5*(p.us1 + p.us2);
		p.rta = 0.5*(p.rt1 + p.rt2);
		p.dea = 0.5*(p.de1 + p.de2);
		p.da  = 0.5*(p.d1  + p.d2 );
		
		if(ityp==3) 
            p.ald = p.dlcon;%//------ increased dissipation length in wake (decrease its reciprocal)
        else
            p.ald = 1.0;
        end
		
		%//----- set and linearize  equilibrium 1/ue due/dx   ...  new  12 oct 94
		if(ityp==2) 
			p.gcc = p.gccon;
			p.hkc     = p.hka - 1.0 - p.gcc/p.rta;
			p.hkc_hka = 1.0;
			p.hkc_rta =  p.gcc/p.rta/p.rta;
			if(p.hkc < 0.01) 
				p.hkc = 0.01;
				p.hkc_hka = 0.0;
				p.hkc_rta = 0.0;
            end
		
        else
			p.gcc = 0.0;
			p.hkc = p.hka - 1.0;
			p.hkc_hka = 1.0;
			p.hkc_rta = 0.0;
        end
		
		p.hr     = p.hkc     / (p.gacon*p.ald*p.hka);
		p.hr_hka = p.hkc_hka / (p.gacon*p.ald*p.hka) - p.hr / p.hka;
		p.hr_rta = p.hkc_rta / (p.gacon*p.ald*p.hka);
		
		p.uq     = (0.5*p.cfa - p.hr*p.hr) / (p.gbcon*p.da);
		p.uq_hka =   -2.0*p.hr*p.hr_hka  / (p.gbcon*p.da);
		p.uq_rta =   -2.0*p.hr*p.hr_rta  / (p.gbcon*p.da);
		p.uq_cfa =   0.5             / (p.gbcon*p.da);
		p.uq_da  = -p.uq/p.da;
		p.uq_upw = p.uq_cfa*(p.cf2-p.cf1) + p.uq_hka*(p.hk2-p.hk1);
		
		p.uq_t1 = (1.0-p.upw)*(p.uq_cfa*p.cf1_t1 + p.uq_hka*p.hk1_t1) + p.uq_upw*p.upw_t1;
		p.uq_d1 = (1.0-p.upw)*(p.uq_cfa*p.cf1_d1 + p.uq_hka*p.hk1_d1) + p.uq_upw*p.upw_d1;
		p.uq_u1 = (1.0-p.upw)*(p.uq_cfa*p.cf1_u1 + p.uq_hka*p.hk1_u1) + p.uq_upw*p.upw_u1;
		p.uq_t2 =       p.upw *(p.uq_cfa*p.cf2_t2 + p.uq_hka*p.hk2_t2) + p.uq_upw*p.upw_t2;
		p.uq_d2 =       p.upw *(p.uq_cfa*p.cf2_d2 + p.uq_hka*p.hk2_d2) + p.uq_upw*p.upw_d2;
		p.uq_u2 =       p.upw *(p.uq_cfa*p.cf2_u2 + p.uq_hka*p.hk2_u2) + p.uq_upw*p.upw_u2;
		p.uq_ms = (1.0-p.upw)*(p.uq_cfa*p.cf1_ms + p.uq_hka*p.hk1_ms) + p.uq_upw*p.upw_ms...
			+       p.upw *(p.uq_cfa*p.cf2_ms + p.uq_hka*p.hk2_ms);
		p.uq_re = (1.0-p.upw)* p.uq_cfa*p.cf1_re + p.upw * p.uq_cfa*p.cf2_re;
		
		p.uq_t1 = p.uq_t1             + 0.5*p.uq_rta*p.rt1_t1;
		p.uq_d1 = p.uq_d1 + 0.5*p.uq_da;
		p.uq_u1 = p.uq_u1             + 0.5*p.uq_rta*p.rt1_u1;
		p.uq_t2 = p.uq_t2             + 0.5*p.uq_rta*p.rt2_t2;
		p.uq_d2 = p.uq_d2 + 0.5*p.uq_da;
		p.uq_u2 = p.uq_u2             + 0.5*p.uq_rta*p.rt2_u2;
		p.uq_ms = p.uq_ms             + 0.5*p.uq_rta*p.rt1_ms...
			+ 0.5*p.uq_rta*p.rt2_ms;
		p.uq_re = p.uq_re             + 0.5*p.uq_rta*p.rt1_re...
			+ 0.5*p.uq_rta*p.rt2_re;
		
		p.scc = p.sccon*1.333/(1.0+p.usa);
		p.scc_usa = -p.scc/(1.0+p.usa);
		
		p.scc_us1 = p.scc_usa*0.5;
		p.scc_us2 = p.scc_usa*0.5;
		
		p.slog = log(p.s2/p.s1);
		p.dxi = p.x2 - p.x1;
		
		p.rezc = p.scc*(p.cqa - p.sa*p.ald)*p.dxi- p.dea*2.0*  p.slog + p.dea*2.0*(p.uq*p.dxi - p.ulog);
		
		p.z_cfa = p.dea*2.0*p.uq_cfa*p.dxi;
		p.z_hka = p.dea*2.0*p.uq_hka*p.dxi;
		p.z_da  = p.dea*2.0*p.uq_da *p.dxi;
		p.z_sl = -p.dea*2.0;
		p.z_ul = -p.dea*2.0;
		p.z_dxi = p.scc    *(p.cqa - p.sa*p.ald)     + p.dea*2.0*p.uq;
		p.z_usa = p.scc_usa*(p.cqa - p.sa*p.ald)*p.dxi;
		p.z_cqa = p.scc*p.dxi;
		p.z_sa = -p.scc*p.dxi*p.ald;
		p.z_dea = 2.0*(p.uq*p.dxi - p.ulog - p.slog);
		
		p.z_upw = p.z_cqa*(p.cq2-p.cq1) + p.z_sa *(p.s2 -p.s1 )+ p.z_cfa*(p.cf2-p.cf1) + p.z_hka*(p.hk2-p.hk1); 
		p.z_de1 = 0.5*p.z_dea;
		p.z_de2 = 0.5*p.z_dea;
		p.z_us1 = 0.5*p.z_usa;
		p.z_us2 = 0.5*p.z_usa;
		p.z_d1  = 0.5*p.z_da;
		p.z_d2  = 0.5*p.z_da;
		p.z_u1  =                 - p.z_ul/p.u1;
		p.z_u2  =                   p.z_ul/p.u2;
		p.z_x1  = -p.z_dxi;
		p.z_x2  =  p.z_dxi;
		p.z_s1  = (1.0-p.upw)*p.z_sa  - p.z_sl/p.s1;
		p.z_s2  =       p.upw *p.z_sa  + p.z_sl/p.s2;
		p.z_cq1 = (1.0-p.upw)*p.z_cqa;
		p.z_cq2 =       p.upw *p.z_cqa;
		p.z_cf1 = (1.0-p.upw)*p.z_cfa;
		p.z_cf2 =       p.upw *p.z_cfa;
		p.z_hk1 = (1.0-p.upw)*p.z_hka;
		p.z_hk2 =       p.upw *p.z_hka;
		
		p.vs1(1,1) = p.z_s1;
		p.vs1(1,2) =        p.z_upw*p.upw_t1 + p.z_de1*p.de1_t1 + p.z_us1*p.us1_t1;
		p.vs1(1,3) = p.z_d1 + p.z_upw*p.upw_d1 + p.z_de1*p.de1_d1 + p.z_us1*p.us1_d1;
		p.vs1(1,4) = p.z_u1 + p.z_upw*p.upw_u1 + p.z_de1*p.de1_u1 + p.z_us1*p.us1_u1;
		p.vs1(1,5) = p.z_x1;
		p.vs2(1,1) = p.z_s2;
		p.vs2(1,2) =        p.z_upw*p.upw_t2 + p.z_de2*p.de2_t2 + p.z_us2*p.us2_t2;
		p.vs2(1,3) = p.z_d2 + p.z_upw*p.upw_d2 + p.z_de2*p.de2_d2 + p.z_us2*p.us2_d2;
		p.vs2(1,4) = p.z_u2 + p.z_upw*p.upw_u2 + p.z_de2*p.de2_u2 + p.z_us2*p.us2_u2;
		p.vs2(1,5) = p.z_x2;
		p.vsm(1)   =        p.z_upw*p.upw_ms + p.z_de1*p.de1_ms + p.z_us1*p.us1_ms+ p.z_de2*p.de2_ms + p.z_us2*p.us2_ms;
		
		p.vs1(1,2) = p.vs1(1,2) + p.z_cq1*p.cq1_t1 + p.z_cf1*p.cf1_t1 + p.z_hk1*p.hk1_t1;
		p.vs1(1,3) = p.vs1(1,3) + p.z_cq1*p.cq1_d1 + p.z_cf1*p.cf1_d1 + p.z_hk1*p.hk1_d1;
		p.vs1(1,4) = p.vs1(1,4) + p.z_cq1*p.cq1_u1 + p.z_cf1*p.cf1_u1 + p.z_hk1*p.hk1_u1;
		
		p.vs2(1,2) = p.vs2(1,2) + p.z_cq2*p.cq2_t2 + p.z_cf2*p.cf2_t2 + p.z_hk2*p.hk2_t2;
		p.vs2(1,3) = p.vs2(1,3) + p.z_cq2*p.cq2_d2 + p.z_cf2*p.cf2_d2 + p.z_hk2*p.hk2_d2;
		p.vs2(1,4) = p.vs2(1,4) + p.z_cq2*p.cq2_u2 + p.z_cf2*p.cf2_u2 + p.z_hk2*p.hk2_u2;
		
		p.vsm(1)   = p.vsm(1)   + p.z_cq1*p.cq1_ms + p.z_cf1*p.cf1_ms + p.z_hk1*p.hk1_ms...
			+ p.z_cq2*p.cq2_ms + p.z_cf2*p.cf2_ms + p.z_hk2*p.hk2_ms;
		p.vsr(1)   =            p.z_cq1*p.cq1_re + p.z_cf1*p.cf1_re...
			+ p.z_cq2*p.cq2_re + p.z_cf2*p.cf2_re;
		p.vsx(1)   = 0.0;
		p.vsrez(1) = -p.rezc;
        end
    end%//endif

	%//**** set up momentum equation
	p.ha = 0.5*(p.h1 + p.h2);
	p.ma = 0.5*(p.m1 + p.m2);
	p.xa = 0.5*(p.x1 + p.x2);
	p.ta = 0.5*(p.t1 + p.t2);
	p.hwa = 0.5*(p.dw1/p.t1 + p.dw2/p.t2);

	%//---- set cf term, using central value cfm for better accuracy in drag
	p.cfx     = 0.50*p.cfm*p.xa/p.ta  +  0.25*(p.cf1*p.x1/p.t1 + p.cf2*p.x2/p.t2);
	p.cfx_xa  = 0.50*p.cfm   /p.ta;
	p.cfx_ta  = -.50*p.cfm*p.xa/p.ta/p.ta;

	p.cfx_x1  = 0.25*p.cf1   /p.t1     + p.cfx_xa*0.5;
	p.cfx_x2  = 0.25*p.cf2   /p.t2     + p.cfx_xa*0.5;
	p.cfx_t1  = -.25*p.cf1*p.x1/p.t1/p.t1  + p.cfx_ta*0.5;
	p.cfx_t2  = -.25*p.cf2*p.x2/p.t2/p.t2  + p.cfx_ta*0.5;
	p.cfx_cf1 = 0.25*    p.x1/p.t1;
	p.cfx_cf2 = 0.25*    p.x2/p.t2;
	p.cfx_cfm = 0.50*    p.xa/p.ta;

	p.btmp = p.ha + 2.0 - p.ma + p.hwa;

	p.rezt  = p.tlog + p.btmp*p.ulog - p.xlog*0.5*p.cfx;
	p.z_cfx = -p.xlog*0.5;
	p.z_ha  =  p.ulog;
	p.z_hwa =  p.ulog;
	p.z_ma  = -p.ulog;
	p.z_xl  =-p.ddlog * 0.5*p.cfx;
	p.z_ul  = p.ddlog * p.btmp;
	p.z_tl  = p.ddlog;

	p.z_cfm = p.z_cfx*p.cfx_cfm;
	p.z_cf1 = p.z_cfx*p.cfx_cf1;
	p.z_cf2 = p.z_cfx*p.cfx_cf2;

	p.z_t1 = -p.z_tl/p.t1 + p.z_cfx*p.cfx_t1 + p.z_hwa*0.5*(-p.dw1/p.t1/p.t1);
	p.z_t2 =  p.z_tl/p.t2 + p.z_cfx*p.cfx_t2 + p.z_hwa*0.5*(-p.dw2/p.t2/p.t2);
	p.z_x1 = -p.z_xl/p.x1 + p.z_cfx*p.cfx_x1;
	p.z_x2 =  p.z_xl/p.x2 + p.z_cfx*p.cfx_x2;
	p.z_u1 = -p.z_ul/p.u1;
	p.z_u2 =  p.z_ul/p.u2;

	p.vs1(2,2) = 0.5*p.z_ha*p.h1_t1 + p.z_cfm*p.cfm_t1 + p.z_cf1*p.cf1_t1 + p.z_t1;
	p.vs1(2,3) = 0.5*p.z_ha*p.h1_d1 + p.z_cfm*p.cfm_d1 + p.z_cf1*p.cf1_d1;
	p.vs1(2,4) = 0.5*p.z_ma*p.m1_u1 + p.z_cfm*p.cfm_u1 + p.z_cf1*p.cf1_u1 + p.z_u1;
	p.vs1(2,5) =                                                p.z_x1;
	p.vs2(2,2) = 0.5*p.z_ha*p.h2_t2 + p.z_cfm*p.cfm_t2 + p.z_cf2*p.cf2_t2 + p.z_t2;
	p.vs2(2,3) = 0.5*p.z_ha*p.h2_d2 + p.z_cfm*p.cfm_d2 + p.z_cf2*p.cf2_d2;
	p.vs2(2,4) = 0.5*p.z_ma*p.m2_u2 + p.z_cfm*p.cfm_u2 + p.z_cf2*p.cf2_u2 + p.z_u2;
	p.vs2(2,5) =                                                p.z_x2;

	p.vsm(2)   = 0.5*p.z_ma*p.m1_ms + p.z_cfm*p.cfm_ms + p.z_cf1*p.cf1_ms...
			 + 0.5*p.z_ma*p.m2_ms + p.z_cf2*p.cf2_ms;
	p.vsr(2)   =                   p.z_cfm*p.cfm_re + p.z_cf1*p.cf1_re...
							   + p.z_cf2*p.cf2_re;
	p.vsx(2)   = 0.0;
	p.vsrez(2) = -p.rezt;

	%//**** set up shape parameter equation

	p.xot1 = p.x1/p.t1;
	p.xot2 = p.x2/p.t2;

	p.ha  = 0.5*(p.h1  + p.h2 );
	p.hsa = 0.5*(p.hs1 + p.hs2);
	p.hca = 0.5*(p.hc1 + p.hc2);
	p.hwa = 0.5*(p.dw1/p.t1 + p.dw2/p.t2);

	p.dix = (1.0-p.upw)*p.di1*p.xot1 + p.upw*p.di2*p.xot2;
	p.cfx = (1.0-p.upw)*p.cf1*p.xot1 + p.upw*p.cf2*p.xot2;
	p.dix_upw = p.di2*p.xot2 - p.di1*p.xot1;
	p.cfx_upw = p.cf2*p.xot2 - p.cf1*p.xot1;

	p.btmp = 2.0*p.hca/p.hsa + 1.0 - p.ha - p.hwa;

	p.rezh  = p.hlog + p.btmp*p.ulog + p.xlog*(0.5*p.cfx-p.dix);
	p.z_cfx =  p.xlog*0.5;
	p.z_dix = -p.xlog;
	p.z_hca = 2.0*p.ulog/p.hsa;
	p.z_ha  = -p.ulog;
	p.z_hwa = -p.ulog;
	p.z_xl  = p.ddlog * (0.5*p.cfx-p.dix);
	p.z_ul  = p.ddlog * p.btmp;
	p.z_hl  = p.ddlog;

	p.z_upw = p.z_cfx*p.cfx_upw + p.z_dix*p.dix_upw;

	p.z_hs1 = -p.hca*p.ulog/p.hsa/p.hsa - p.z_hl/p.hs1;
	p.z_hs2 = -p.hca*p.ulog/p.hsa/p.hsa + p.z_hl/p.hs2;

	p.z_cf1 = (1.0-p.upw)*p.z_cfx*p.xot1;
	p.z_cf2 =      p.upw *p.z_cfx*p.xot2;
	p.z_di1 = (1.0-p.upw)*p.z_dix*p.xot1;
	p.z_di2 =      p.upw *p.z_dix*p.xot2;

	p.z_t1 = (1.0-p.upw)*(p.z_cfx*p.cf1 + p.z_dix*p.di1)*(-p.xot1/p.t1);
	p.z_t2 =      p.upw *(p.z_cfx*p.cf2 + p.z_dix*p.di2)*(-p.xot2/p.t2);
	p.z_x1 = (1.0-p.upw)*(p.z_cfx*p.cf1 + p.z_dix*p.di1)/ p.t1        - p.z_xl/p.x1;
	p.z_x2 =      p.upw *(p.z_cfx*p.cf2 + p.z_dix*p.di2)/ p.t2        + p.z_xl/p.x2;
	p.z_u1 =                                              - p.z_ul/p.u1;
	p.z_u2 =                                                p.z_ul/p.u2;

	p.z_t1 = p.z_t1 + p.z_hwa*0.5*(-p.dw1/p.t1/p.t1);
	p.z_t2 = p.z_t2 + p.z_hwa*0.5*(-p.dw2/p.t2/p.t2);

	p.vs1(3,1) =                               p.z_di1*p.di1_s1;
	p.vs1(3,2) = p.z_hs1*p.hs1_t1 + p.z_cf1*p.cf1_t1 + p.z_di1*p.di1_t1 + p.z_t1;
	p.vs1(3,3) = p.z_hs1*p.hs1_d1 + p.z_cf1*p.cf1_d1 + p.z_di1*p.di1_d1;
	p.vs1(3,4) = p.z_hs1*p.hs1_u1 + p.z_cf1*p.cf1_u1 + p.z_di1*p.di1_u1 + p.z_u1;
	p.vs1(3,5) =                                              p.z_x1;
	p.vs2(3,1) =                               p.z_di2*p.di2_s2;
	p.vs2(3,2) = p.z_hs2*p.hs2_t2 + p.z_cf2*p.cf2_t2 +p.z_di2*p.di2_t2 + p.z_t2;
	p.vs2(3,3) = p.z_hs2*p.hs2_d2 + p.z_cf2*p.cf2_d2 + p.z_di2*p.di2_d2;
	p.vs2(3,4) = p.z_hs2*p.hs2_u2 + p.z_cf2*p.cf2_u2 + p.z_di2*p.di2_u2 + p.z_u2;
	p.vs2(3,5) =                                              p.z_x2;
	p.vsm(3)   = p.z_hs1*p.hs1_ms + p.z_cf1*p.cf1_ms + p.z_di1*p.di1_ms...
	+ p.z_hs2*p.hs2_ms + p.z_cf2*p.cf2_ms + p.z_di2*p.di2_ms;
	p.vsr(3)   = p.z_hs1*p.hs1_re + p.z_cf1*p.cf1_re + p.z_di1*p.di1_re...
	+ p.z_hs2*p.hs2_re + p.z_cf2*p.cf2_re + p.z_di2*p.di2_re;

	p.vs1(3,2) = p.vs1(3,2) + 0.5*(p.z_hca*p.hc1_t1+p.z_ha*p.h1_t1) + p.z_upw*p.upw_t1;
	p.vs1(3,3) = p.vs1(3,3) + 0.5*(p.z_hca*p.hc1_d1+p.z_ha*p.h1_d1) + p.z_upw*p.upw_d1;
	p.vs1(3,4) = p.vs1(3,4) + 0.5*(p.z_hca*p.hc1_u1           ) + p.z_upw*p.upw_u1;
	p.vs2(3,2) = p.vs2(3,2) + 0.5*(p.z_hca*p.hc2_t2+p.z_ha*p.h2_t2) + p.z_upw*p.upw_t2;
	p.vs2(3,3) = p.vs2(3,3) + 0.5*(p.z_hca*p.hc2_d2+p.z_ha*p.h2_d2) + p.z_upw*p.upw_d2;
	p.vs2(3,4) = p.vs2(3,4) + 0.5*(p.z_hca*p.hc2_u2           ) + p.z_upw*p.upw_u2;

	p.vsm(3)   = p.vsm(3)   + 0.5*(p.z_hca*p.hc1_ms           ) + p.z_upw*p.upw_ms...
	+ 0.5*(p.z_hca*p.hc2_ms           );

	p.vsx(3)   = 0.0;
	p.vsrez(3) = -p.rezh;

% 	return true;
% }
end
function [p]=blmid(ityp,p)
% ////--------------------------------------
% //     blmid function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% 	//---------------------------------------------------- 
% 	//     calculates midpoint skin friction p.cfm
% 	//
% 	//      ityp = 1 :  laminar
% 	//      ityp = 2 :  turbulent
% 	//      ityp = 3 :  turbulent wake
% 	//----------------------------------------------------
% 	//
% 	double p.hka, p.rta, p.ma, p.cfm_rta, p.cfm_ma;
% 	double p.cfml, p.cfml_hka, p.cfml_rta, p.cfml_ma, p.cfm_hka;
% 	//---- set similarity variables if not defined
	if(p.simi) 
		p.hk1    = p.hk2;
		p.hk1_t1 = p.hk2_t2;
		p.hk1_d1 = p.hk2_d2;
		p.hk1_u1 = p.hk2_u2;
		p.hk1_ms = p.hk2_ms;
		p.rt1    = p.rt2;
		p.rt1_t1 = p.rt2_t2;
		p.rt1_u1 = p.rt2_u2;
		p.rt1_ms = p.rt2_ms;
		p.rt1_re = p.rt2_re;
		p.m1    = p.m2;
		p.m1_u1 = p.m2_u2;
		p.m1_ms = p.m2_ms;
    end
	
% 	//---- define stuff for midpoint cf
	p.hka = 0.5*(p.hk1 + p.hk2);
	p.rta = 0.5*(p.rt1 + p.rt2);
	p.ma  = 0.5*(p.m1  + p.m2 );
	
% 	//---- midpoint skin friction coefficient  (zero in wake)
	if(ityp==3)
		p.cfm     = 0.0;
		p.cfm_hka = 0.0;
		p.cfm_rta = 0.0;
		p.cfm_ma  = 0.0;
		p.cfm_ms  = 0.0;
	
	else 
		if(ityp==1) 
            [p.cfm, p.cfm_hka, p.cfm_rta, p.cfm_ma]=cfl( p.hka, p.rta );
		else 
			[p.cfm, p.cfm_hka, p.cfm_rta, p.cfm_ma]=cft( p.hka, p.rta, p.ma );
			[p.cfml, p.cfml_hka, p.cfml_rta, p.cfml_ma]=cfl( p.hka, p.rta);
			if(p.cfml>p.cfm) 
% 				//ccc      write(*,*) 'cft cfl rt hk:', p.cfm, p.cfml, p.rta, p.hka, 0.5*(x1+x2) 
				p.cfm     = p.cfml;
				p.cfm_hka = p.cfml_hka;
				p.cfm_rta = p.cfml_rta;
				p.cfm_ma  = p.cfml_ma;
            end
        end
    end
	p.cfm_u1 = 0.5*(p.cfm_hka*p.hk1_u1 + p.cfm_ma*p.m1_u1 + p.cfm_rta*p.rt1_u1);
	p.cfm_t1 = 0.5*(p.cfm_hka*p.hk1_t1 +                p.cfm_rta*p.rt1_t1);
	p.cfm_d1 = 0.5*(p.cfm_hka*p.hk1_d1                                );
	
	p.cfm_u2 = 0.5*(p.cfm_hka*p.hk2_u2 + p.cfm_ma*p.m2_u2 + p.cfm_rta*p.rt2_u2);
	p.cfm_t2 = 0.5*(p.cfm_hka*p.hk2_t2 +                p.cfm_rta*p.rt2_t2);
	p.cfm_d2 = 0.5*(p.cfm_hka*p.hk2_d2                                );
	
	p.cfm_ms = 0.5*(p.cfm_hka*p.hk1_ms + p.cfm_ma*p.m1_ms + p.cfm_rta*p.rt1_ms...
		+ p.cfm_hka*p.hk2_ms + p.cfm_ma*p.m2_ms + p.cfm_rta*p.rt2_ms);
	p.cfm_re = 0.5*(                                p.cfm_rta*p.rt1_re...
		+ p.cfm_rta*p.rt2_re);
	
% 	return true;
% }
end
function [p]=blvar(ityp,p)
% ////--------------------------------------
% //     blvar function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% //----------------------------------------------------
% //     calculates all secondary "2" variables from
% //     the primary "2" variables x2, u2, p.t2, p.d2, p.s2.
% //     also calculates the sensitivities of the
% //     secondary variables wrt the primary variables.
% //
% % //      ityp = 1 :  laminar
% //      ityp = 2 :  turbulent
% //      ityp = 3 :  turbulent wake
% //----------------------------------------------------
% 
% 	double p.hs2_hk2, p.hs2_rt2, p.hs2_m2;
% 	double p.hc2_hk2, p.hc2_m2, p.us2_hs2, p.us2_hk2, p.us2_h2;
% 	double p.gcc, p.hkc, p.hkc_hk2, p.hkc_rt2, p.hkb, p.usb;
% 	double p.cq2_hs2, p.cq2_us2, p.cq2_hk2;
% 	double p.cq2_rt2, p.cq2_h2, p.cf2_hk2, p.cf2_rt2, p.cf2_m2;
% 	double p.cf2l, p.cf2l_hk2, p.cf2l_rt2, p.cf2l_m2;
% 	double p.di2_hk2, p.di2_rt2, p.cf2t, p.cf2t_hk2;
% 	double p.cf2t_rt2, p.cf2t_m2, p.cf2t_u2, p.cf2t_t2;
% 	double p.cf2t_d2, p.cf2t_ms, p.cf2t_re, p.di2_hs2;
% 	double p.di2_us2, p.di2_cf2t, p.hmin, p.hm_rt2;
% 	double p.grt, p.fl, p.fl_hk2, p.fl_rt2, p.tfl;
% 	double p.dfac, p.df_fl, p.df_hk2, p.df_rt2, p.dd, p.dd_hs2;
% 	double p.dd_us2, p.dd_s2, p.dd_rt2, p.di2l, p.di2l_hk2, p.di2l_rt2, p.de2_hk2, p.hdmax;
% 	
% 	
% //	double p.gbcon, p.gccon, p.ctcon, hkc2;//were are they initialized ?

    
	if(ityp==3) 
        p.hk2 = max(p.hk2,1.00005);
    end
	if(ityp~=3) 
        p.hk2 = max(p.hk2,1.05000);
    end
% 	//---- density thickness shape parameter     ( h** )
% [hc,hc_hk, hc_msq]=hct( hk,  msq)
	[p.hc2, p.hc2_hk2, p.hc2_m2]=hct( p.hk2, p.m2 );
	p.hc2_u2 = p.hc2_hk2*p.hk2_u2 + p.hc2_m2*p.m2_u2;
	p.hc2_t2 = p.hc2_hk2*p.hk2_t2;
	p.hc2_d2 = p.hc2_hk2*p.hk2_d2;
	p.hc2_ms = p.hc2_hk2*p.hk2_ms + p.hc2_m2*p.m2_ms;
	
% 	//---- set ke thickness shape parameter from  h - h*  correlations
	if(ityp==1)
        [p.hs2,p.hs2_hk2, p.hs2_rt2, p.hs2_m2]=hsl(p.hk2 );
    else
        [p.hs2, p.hs2_hk2, p.hs2_rt2, p.hs2_m2]=hst(p.hk2, p.rt2, p.m2 );
    end
	
	
	p.hs2_u2 = p.hs2_hk2*p.hk2_u2 + p.hs2_rt2*p.rt2_u2 + p.hs2_m2*p.m2_u2;
	p.hs2_t2 = p.hs2_hk2*p.hk2_t2 + p.hs2_rt2*p.rt2_t2;
	p.hs2_d2 = p.hs2_hk2*p.hk2_d2;
	p.hs2_ms = p.hs2_hk2*p.hk2_ms + p.hs2_rt2*p.rt2_ms + p.hs2_m2*p.m2_ms;
	p.hs2_re =                  p.hs2_rt2*p.rt2_re;
	
% 	//---- normalized slip velocity  us
	p.us2     = 0.5*p.hs2*( 1.0 - (p.hk2-1.0)/(p.gbcon*p.h2) );
	p.us2_hs2 = 0.5  *  ( 1.0 - (p.hk2-1.0)/(p.gbcon*p.h2) );
	p.us2_hk2 = 0.5*p.hs2*(      -  1.0     /(p.gbcon*p.h2) );
	p.us2_h2  = 0.5*p.hs2*         (p.hk2-1.0)/(p.gbcon*p.h2*p.h2);
	
	p.us2_u2 = p.us2_hs2*p.hs2_u2 + p.us2_hk2*p.hk2_u2;
	p.us2_t2 = p.us2_hs2*p.hs2_t2 + p.us2_hk2*p.hk2_t2 + p.us2_h2*p.h2_t2;
	p.us2_d2 = p.us2_hs2*p.hs2_d2 + p.us2_hk2*p.hk2_d2 + p.us2_h2*p.h2_d2;
	p.us2_ms = p.us2_hs2*p.hs2_ms + p.us2_hk2*p.hk2_ms;
	p.us2_re = p.us2_hs2*p.hs2_re;
	
	if(ityp<=2 && p.us2>0.95) 
% 		//       write(*,*) 'blvar: us clamped:', p.us2
		p.us2 = 0.98;
		p.us2_u2 = 0.0;
		p.us2_t2 = 0.0;
		p.us2_d2 = 0.0;
		p.us2_ms = 0.0;
		p.us2_re = 0.0;
    end
	
	if(ityp==3 && p.us2>0.99995) 
% 		//       write(*,*) 'blvar: wake us clamped:', p.us2
		p.us2 = 0.99995;
		p.us2_u2 = 0.0;
		p.us2_t2 = 0.0;
		p.us2_d2 = 0.0;
		p.us2_ms = 0.0;
		p.us2_re = 0.0;
    end
	
% 	//---- equilibrium wake layer shear coefficient (ctau)eq ** 1/2
% 	//   ...  new  12 oct 94
	p.gcc = 0.0;
	p.hkc = p.hk2 - 1.0;
	p.hkc_hk2 = 1.0;
	p.hkc_rt2 = 0.0;
	if(ityp==2) 
		p.gcc = p.gccon;
		p.hkc     = p.hk2 - 1.0 - p.gcc/p.rt2;
		p.hkc_hk2 = 1.0;
		p.hkc_rt2 =             p.gcc/p.rt2/p.rt2;
		if(p.hkc < 0.01) 
			p.hkc = 0.01;
			p.hkc_hk2 = 0.0;
			p.hkc_rt2 = 0.0;
        end
    end
	
	p.hkb = p.hk2 - 1.0;
	p.usb = 1.0 - p.us2;
	p.cq2     =...
		sqrt( p.ctcon*p.hs2*p.hkb*p.hkc*p.hkc / (p.usb*p.h2*p.hk2*p.hk2) );
	p.cq2_hs2 = p.ctcon    *p.hkb*p.hkc*p.hkc / (p.usb*p.h2*p.hk2*p.hk2)       * 0.5/p.cq2;
	p.cq2_us2 = p.ctcon*p.hs2*p.hkb*p.hkc*p.hkc / (p.usb*p.h2*p.hk2*p.hk2) / p.usb * 0.5/p.cq2;
	p.cq2_hk2 = p.ctcon*p.hs2    *p.hkc*p.hkc / (p.usb*p.h2*p.hk2*p.hk2)       * 0.5/p.cq2...
	- p.ctcon*p.hs2*p.hkb*p.hkc*p.hkc / (p.usb*p.h2*p.hk2*p.hk2*p.hk2) * 2.0 * 0.5/p.cq2...
		+ p.ctcon*p.hs2*p.hkb*p.hkc     / (p.usb*p.h2*p.hk2*p.hk2) * 2.0 * 0.5/p.cq2...
		*p.hkc_hk2;
	p.cq2_rt2 = p.ctcon*p.hs2*p.hkb*p.hkc    / (p.usb*p.h2*p.hk2*p.hk2) * 2.0 * 0.5/p.cq2...
		*p.hkc_rt2;
	p.cq2_h2  =-p.ctcon*p.hs2*p.hkb*p.hkc*p.hkc / (p.usb*p.h2*p.hk2*p.hk2) / p.h2  * 0.5/p.cq2;
	
	p.cq2_u2 = p.cq2_hs2*p.hs2_u2 + p.cq2_us2*p.us2_u2 + p.cq2_hk2*p.hk2_u2;
	p.cq2_t2 = p.cq2_hs2*p.hs2_t2 + p.cq2_us2*p.us2_t2 + p.cq2_hk2*p.hk2_t2;
	p.cq2_d2 = p.cq2_hs2*p.hs2_d2 + p.cq2_us2*p.us2_d2 + p.cq2_hk2*p.hk2_d2;
	p.cq2_ms = p.cq2_hs2*p.hs2_ms + p.cq2_us2*p.us2_ms + p.cq2_hk2*p.hk2_ms;
	p.cq2_re = p.cq2_hs2*p.hs2_re + p.cq2_us2*p.us2_re;
	
	p.cq2_u2 = p.cq2_u2                + p.cq2_rt2*p.rt2_u2;
	p.cq2_t2 = p.cq2_t2 + p.cq2_h2*p.h2_t2 + p.cq2_rt2*p.rt2_t2;
	p.cq2_d2 = p.cq2_d2 + p.cq2_h2*p.h2_d2;
	p.cq2_ms = p.cq2_ms                + p.cq2_rt2*p.rt2_ms;
	p.cq2_re = p.cq2_re                + p.cq2_rt2*p.rt2_re;
	
% 	//---- set skin friction coefficient 
	if(ityp==3) 
% 		//----- wake
		p.cf2     = 0.0;
		p.cf2_hk2 = 0.0;
		p.cf2_rt2 = 0.0;
		p.cf2_m2  = 0.0;
	
    else
		if(ityp==1) 
% 			//----- laminar
			[ p.cf2, p.cf2_hk2, p.cf2_rt2, p.cf2_m2]=cfl(p.hk2, p.rt2);
        else
% 			//----- turbulent
			[p.cf2, p.cf2_hk2, p.cf2_rt2, p.cf2_m2]=cft(p.hk2, p.rt2, p.m2);
			[p.cf2l, p.cf2l_hk2, p.cf2l_rt2, p.cf2l_m2]=cfl(p.hk2, p.rt2);
			if(p.cf2l>p.cf2) 
% 				//------- laminar cf is greater than turbulent cf -- use laminar
% 				//-       (this will only occur for unreasonably small rtheta)
				p.cf2     = p.cf2l;
				p.cf2_hk2 = p.cf2l_hk2;
				p.cf2_rt2 = p.cf2l_rt2;
				p.cf2_m2  = p.cf2l_m2;
            end
        end
    end
	
	p.cf2_u2 = p.cf2_hk2*p.hk2_u2 + p.cf2_rt2*p.rt2_u2 + p.cf2_m2*p.m2_u2;
	p.cf2_t2 = p.cf2_hk2*p.hk2_t2 + p.cf2_rt2*p.rt2_t2;
	p.cf2_d2 = p.cf2_hk2*p.hk2_d2;
	p.cf2_ms = p.cf2_hk2*p.hk2_ms + p.cf2_rt2*p.rt2_ms + p.cf2_m2*p.m2_ms;
	p.cf2_re =                  p.cf2_rt2*p.rt2_re;
	
% 	//---- dissipation function    2 cd / h*
	if(ityp==1) 
		
% 		//----- laminar
		[p.di2, p.di2_hk2, p.di2_rt2]=dil( p.hk2, p.rt2 );
		
		p.di2_u2 = p.di2_hk2*p.hk2_u2 + p.di2_rt2*p.rt2_u2;
		p.di2_t2 = p.di2_hk2*p.hk2_t2 + p.di2_rt2*p.rt2_t2;
		p.di2_d2 = p.di2_hk2*p.hk2_d2;
		p.di2_s2 = 0.0;
		p.di2_ms = p.di2_hk2*p.hk2_ms + p.di2_rt2*p.rt2_ms;
		p.di2_re =                  p.di2_rt2*p.rt2_re;
	
    else
		if(ityp==2) 
			
			
% 			//----- turbulent wall contribution
			[p.cf2t, p.cf2t_hk2, p.cf2t_rt2, p.cf2t_m2]=cft(p.hk2, p.rt2, p.m2);
			p.cf2t_u2 = p.cf2t_hk2*p.hk2_u2 + p.cf2t_rt2*p.rt2_u2 + p.cf2t_m2*p.m2_u2;
			p.cf2t_t2 = p.cf2t_hk2*p.hk2_t2 + p.cf2t_rt2*p.rt2_t2;
			p.cf2t_d2 = p.cf2t_hk2*p.hk2_d2;
			p.cf2t_ms = p.cf2t_hk2*p.hk2_ms + p.cf2t_rt2*p.rt2_ms + p.cf2t_m2*p.m2_ms;
			p.cf2t_re =                   p.cf2t_rt2*p.rt2_re;
			
			p.di2      =  ( 0.5*p.cf2t*p.us2 ) * 2.0/p.hs2;
			p.di2_hs2  = -( 0.5*p.cf2t*p.us2 ) * 2.0/p.hs2/p.hs2;
			p.di2_us2  =  ( 0.5*p.cf2t     ) * 2.0/p.hs2;
			p.di2_cf2t =  ( 0.5     *p.us2 ) * 2.0/p.hs2;
			
			p.di2_s2 = 0.0;
			p.di2_u2 = p.di2_hs2*p.hs2_u2 + p.di2_us2*p.us2_u2 + p.di2_cf2t*p.cf2t_u2;
			p.di2_t2 = p.di2_hs2*p.hs2_t2 + p.di2_us2*p.us2_t2 + p.di2_cf2t*p.cf2t_t2;
			p.di2_d2 = p.di2_hs2*p.hs2_d2 + p.di2_us2*p.us2_d2 + p.di2_cf2t*p.cf2t_d2;
			p.di2_ms = p.di2_hs2*p.hs2_ms + p.di2_us2*p.us2_ms + p.di2_cf2t*p.cf2t_ms;
			p.di2_re = p.di2_hs2*p.hs2_re + p.di2_us2*p.us2_re + p.di2_cf2t*p.cf2t_re;
			
% 			//----- set minimum hk for wake layer to still exist
			p.grt = log(p.rt2);
			p.hmin = 1.0 + 2.1/p.grt;
			p.hm_rt2 = -(2.1/p.grt/p.grt) / p.rt2;
			
% 			//----- set factor p.dfac for correcting wall dissipation for very low hk
			p.fl = (p.hk2-1.0)/(p.hmin-1.0);
			p.fl_hk2 =   1.0/(p.hmin-1.0);
			p.fl_rt2 = ( -p.fl/(p.hmin-1.0) ) * p.hm_rt2;
			
			p.tfl = tanh(p.fl);
			p.dfac  = 0.5 + 0.5* p.tfl;
			p.df_fl =       0.5*(1.0 - p.tfl*p.tfl);
			
			p.df_hk2 = p.df_fl*p.fl_hk2;
			p.df_rt2 = p.df_fl*p.fl_rt2;
			
			p.di2_s2 = p.di2_s2*p.dfac;
			p.di2_u2 = p.di2_u2*p.dfac + p.di2*(p.df_hk2*p.hk2_u2 + p.df_rt2*p.rt2_u2);
			p.di2_t2 = p.di2_t2*p.dfac + p.di2*(p.df_hk2*p.hk2_t2 + p.df_rt2*p.rt2_t2);
			p.di2_d2 = p.di2_d2*p.dfac + p.di2*(p.df_hk2*p.hk2_d2                );
			p.di2_ms = p.di2_ms*p.dfac + p.di2*(p.df_hk2*p.hk2_ms + p.df_rt2*p.rt2_ms);
			p.di2_re = p.di2_re*p.dfac + p.di2*(                p.df_rt2*p.rt2_re);
			p.di2    = p.di2   *p.dfac;
		
        else
			
% 			//----- zero wall contribution for wake
			p.di2    = 0.0;
			p.di2_s2 = 0.0;
			p.di2_u2 = 0.0;
			p.di2_t2 = 0.0;
			p.di2_d2 = 0.0;
			p.di2_ms = 0.0;
			p.di2_re = 0.0;
			
        end
    end	
% 	//---- add on turbulent outer layer contribution
	if(ityp~=1) 
		
		p.dd     =  p.s2*p.s2 *  (0.995-p.us2) * 2.0/p.hs2;
		p.dd_hs2 = -p.s2*p.s2 *  (0.995-p.us2) * 2.0/p.hs2/p.hs2;
		p.dd_us2 = -p.s2*p.s2               * 2.0/p.hs2;
		p.dd_s2  =  p.s2*2.0* (0.995-p.us2) * 2.0/p.hs2;
		
		p.di2    = p.di2    + p.dd;
		p.di2_s2 =          p.dd_s2;
		p.di2_u2 = p.di2_u2 + p.dd_hs2*p.hs2_u2 + p.dd_us2*p.us2_u2;
		p.di2_t2 = p.di2_t2 + p.dd_hs2*p.hs2_t2 + p.dd_us2*p.us2_t2;
		p.di2_d2 = p.di2_d2 + p.dd_hs2*p.hs2_d2 + p.dd_us2*p.us2_d2;
		p.di2_ms = p.di2_ms + p.dd_hs2*p.hs2_ms + p.dd_us2*p.us2_ms;
		p.di2_re = p.di2_re + p.dd_hs2*p.hs2_re + p.dd_us2*p.us2_re;
		
% 		//----- add laminar stress contribution to outer layer cd
		p.dd     =  0.15*(0.995-p.us2)*(0.995-p.us2) / p.rt2  * 2.0/p.hs2;
		p.dd_us2 = -0.15*(0.995-p.us2)*2.0 / p.rt2  * 2.0/p.hs2;
		p.dd_hs2 = -p.dd/p.hs2;
		p.dd_rt2 = -p.dd/p.rt2;
		
		p.di2    = p.di2    + p.dd;
		p.di2_u2 = p.di2_u2 + p.dd_hs2*p.hs2_u2 + p.dd_us2*p.us2_u2 + p.dd_rt2*p.rt2_u2;
		p.di2_t2 = p.di2_t2 + p.dd_hs2*p.hs2_t2 + p.dd_us2*p.us2_t2 + p.dd_rt2*p.rt2_t2;
		p.di2_d2 = p.di2_d2 + p.dd_hs2*p.hs2_d2 + p.dd_us2*p.us2_d2;
		p.di2_ms = p.di2_ms + p.dd_hs2*p.hs2_ms + p.dd_us2*p.us2_ms + p.dd_rt2*p.rt2_ms;
		p.di2_re = p.di2_re + p.dd_hs2*p.hs2_re + p.dd_us2*p.us2_re + p.dd_rt2*p.rt2_re;
		
    end
	
	if(ityp==2) 
		[p.di2l, p.di2l_hk2, p.di2l_rt2]=dil( p.hk2, p.rt2 );
		
		if(p.di2l>p.di2) 
% 			//------- laminar cd is greater than turbulent cd -- use laminar
% 			//-       (this will only occur for unreasonably small rtheta)
			p.di2    = p.di2l;
			p.di2_s2 = 0.0;
			p.di2_u2 = p.di2l_hk2*p.hk2_u2 + p.di2l_rt2*p.rt2_u2;
			p.di2_t2 = p.di2l_hk2*p.hk2_t2 + p.di2l_rt2*p.rt2_t2;
			p.di2_d2 = p.di2l_hk2*p.hk2_d2;
			p.di2_ms = p.di2l_hk2*p.hk2_ms + p.di2l_rt2*p.rt2_ms;
			p.di2_re =                   p.di2l_rt2*p.rt2_re;
        end
    end
	
	if(ityp==3) 
% 		//------ laminar wake cd
		[p.di2l, p.di2l_hk2, p.di2l_rt2]=dilw( p.hk2, p.rt2 );
		if(p.di2l > p.di2) 
% 			//------- laminar wake cd is greater than turbulent cd -- use laminar
% 			//-       (this will only occur for unreasonably small rtheta)
			p.di2    = p.di2l;
			p.di2_s2 = 0.0;
			p.di2_u2 = p.di2l_hk2*p.hk2_u2 + p.di2l_rt2*p.rt2_u2;
			p.di2_t2 = p.di2l_hk2*p.hk2_t2 + p.di2l_rt2*p.rt2_t2;
			p.di2_d2 = p.di2l_hk2*p.hk2_d2;
			p.di2_ms = p.di2l_hk2*p.hk2_ms + p.di2l_rt2*p.rt2_ms;
			p.di2_re =                   p.di2l_rt2*p.rt2_re;
        end
    end
	
	if(ityp==3) 
% 		//----- double dissipation for the wake (two wake halves)
		p.di2    = p.di2   *2.0;
		p.di2_s2 = p.di2_s2*2.0;
		p.di2_u2 = p.di2_u2*2.0;
		p.di2_t2 = p.di2_t2*2.0;
		p.di2_d2 = p.di2_d2*2.0;
		p.di2_ms = p.di2_ms*2.0;
		p.di2_re = p.di2_re*2.0;
    end
	
% 	//---- bl thickness (delta) from simplified green's correlation
	p.de2     = (3.15 + 1.72/(p.hk2-1.0)   )*p.t2  +  p.d2;
	p.de2_hk2 = (     - 1.72/(p.hk2-1.0)/(p.hk2-1.0))*p.t2;
	
	p.de2_u2 = p.de2_hk2*p.hk2_u2;
	p.de2_t2 = p.de2_hk2*p.hk2_t2 + (3.15 + 1.72/(p.hk2-1.0));
	p.de2_d2 = p.de2_hk2*p.hk2_d2 + 1.0;
	p.de2_ms = p.de2_hk2*p.hk2_ms;
	
	
	p.hdmax = 12.0;
	if(p.de2 > p.hdmax*p.t2) 
		
		p.de2    = p.hdmax*p.t2;
		p.de2_u2 =  0.0;
		p.de2_t2 = p.hdmax;
		p.de2_d2 =  0.0;
		p.de2_ms =  0.0;
    end
	
% 	return true;
% }
end
function [hc,hc_hk, hc_msq]=hct( hk,  msq)
% ////--------------------------------------
% //     hct function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% {
% //---- density shape parameter    (from whitfield)
	hc     = msq * (0.064/(hk-0.8) + 0.251);
	hc_hk  = msq * (-.064/(hk-0.8)/(hk-0.8));
	hc_msq =        0.064/(hk-0.8) + 0.251;
	
% 	return true;
% }
end
function [hs,hs_hk, hs_rt, hs_msq]=hsl(hk)
% ////--------------------------------------
% //     hsl function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% //---- laminar hs correlation
	if(hk<4.35) 
		tmp = hk - 4.35;
		hs    = 0.0111*tmp*tmp/(hk+1.0)...
			- 0.0278*tmp*tmp*tmp/(hk+1.0)  + 1.528...
			- 0.0002*(tmp*hk)*(tmp*hk);
		hs_hk = 0.0111*(2.0*tmp    - tmp*tmp/(hk+1.0))/(hk+1.0)...
			- 0.0278*(3.0*tmp*tmp - tmp*tmp*tmp/(hk+1.0))/(hk+1.0)...
			- 0.0002*2.0*tmp*hk * (tmp + hk);
	
    else
		hs    = 0.015*    (hk-4.35)*(hk-4.35)/hk + 1.528;
		hs_hk = 0.015*2.0*(hk-4.35)   /hk...
		- 0.015*    (hk-4.35)* (hk-4.35)/hk/hk;
    end
	
	hs_rt  = 0.0;
	hs_msq = 0.0;
	
% 	return true;
% }
end
function [hs, hs_hk, hs_rt, hs_msq]=hst(hk, rt, msq)
% ////--------------------------------------
% //     hst function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% //---- turbulent hs correlation
     
	 hsmin = 1.5;
	 dhsinf = 0.015;
% 	double rtz, rtz_rt, ho, ho_rt;

% //---- ###  12/4/94
% //---- limited rtheta dependence for rtheta < 200

	if(rt>400.0) 
		ho    = 3.0 + 400.0/rt;
		ho_rt =      - 400.0/rt/rt;
	
    else
		ho    = 4.0;
		ho_rt = 0.0;
    end
	
	if(rt>200.0) 
		rtz    = rt;
		rtz_rt = 1.0;
	
    else
		rtz    = 200.0;
		rtz_rt = 0.0;
    end

	if(hk<ho) 
% 		//----- attached branch
% 		//----- new correlation  29 nov 91
% 		//-     (from  arctan(y+) + schlichting  profiles)
		hr    =   (ho - hk)/(ho-1.0);
		 hr_hk =      - 1.0/(ho-1.0);
		 hr_rt = (1.0 - hr)/(ho-1.0) * ho_rt;
		hs    = (2.0-hsmin-4.0/rtz)*hr*hr  * 1.5/(hk+0.5) + hsmin...
			+ 4.0/rtz;
		hs_hk =-(2.0-hsmin-4.0/rtz)*hr*hr  * 1.5/(hk+0.5)/(hk+0.5)...
			+ (2.0-hsmin-4.0/rtz)*hr*2.0 * 1.5/(hk+0.5) * hr_hk;
		hs_rt = (2.0-hsmin-4.0/rtz)*hr*2.0 * 1.5/(hk+0.5) * hr_rt...
			+ (hr*hr * 1.5/(hk+0.5) - 1.0)*4.0/rtz/rtz * rtz_rt;
	
    else
		
% 		//----- separated branch
		 grt = log(rtz);
		 hdif = hk - ho ;
		 rtmp = hk - ho + 4.0/grt;
		 htmp    = 0.007*grt/rtmp/rtmp + dhsinf/hk;
		 htmp_hk = -.014*grt/rtmp/rtmp/rtmp - dhsinf/hk/hk;
		 htmp_rt = -.014*grt/rtmp/rtmp/rtmp * (-ho_rt - 4.0/grt/grt/rtz * rtz_rt)...
						+ 0.007  /rtmp/rtmp / rtz * rtz_rt;
		hs    = hdif*hdif * htmp + hsmin + 4.0/rtz;
		hs_hk = hdif*2.0* htmp + hdif*hdif * htmp_hk;
		hs_rt = hdif*hdif * htmp_rt      - 4.0/rtz/rtz * rtz_rt...
			    + hdif*2.0* htmp * (-ho_rt);
		
    end

% //---- whitfield's minor additional compressibility correction
      fm = 1.0 + 0.014*msq;
      hs     = ( hs + 0.028*msq) / fm;
      hs_hk  = ( hs_hk          ) / fm;
      hs_rt  = ( hs_rt          ) / fm;
      hs_msq = 0.028/fm - 0.014*(hs)/fm;

%       return true;
      
% }
end
function [cf, cf_hk, cf_rt, cf_msq]=cfl(hk, rt)
% ////--------------------------------------
% //     cfl function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% 				double &cf, double &cf_hk, double &cf_rt, double &cf_msq){
% 	//---- laminar skin friction function  ( cf )    ( from falkner-skan )
% 	double tmp;
	if(hk<5.5) 
		tmp = (5.5-hk)*(5.5-hk)*(5.5-hk) / (hk+1.0);
		cf    = ( 0.0727*tmp                      - 0.07       )/rt;
		cf_hk = ( -.0727*tmp*3.0/(5.5-hk) - 0.0727*tmp/(hk+1.0))/rt;
	
    else
		tmp = 1.0 - 1.0/(hk-4.5);
		cf    = ( 0.015*tmp*tmp      - 0.07  ) / rt;
		cf_hk = ( 0.015*tmp*2.0/(hk-4.5)/(hk-4.5) ) / rt;
    end
	cf_rt = -cf/rt;
	cf_msq = 0.0;
	
% 	return true;
% }
end
function [cf, cf_hk, cf_rt, cf_msq]=cft( hk,  rt,  msq) 
% ////--------------------------------------
% //     cft function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

	gam =1.4;
% 	double gm1, f_arg, fc, grt, gex, thk, cfo;
	
% 	//---- turbulent skin friction function  ( cf )    (coles)
	gm1 = gam - 1.0;
	fc  = sqrt(1.0 + 0.5*gm1*msq);
	grt = log(rt/fc);
	grt = max(grt,3.0);
	
	gex = -1.74 - 0.31*hk;
	
	f_arg = -1.33*hk;
	f_arg = max(-20.0, f_arg );
	
	thk = tanh(4.0 - hk/0.875);
	
	cfo =  0.3*exp(f_arg) * ((grt/2.3026)^gex);
	cf     = ( cfo  +  0.00011*(thk-1.0) ) / fc;
	cf_hk  = (-1.33*cfo - 0.31*log(grt/2.3026)*cfo...
		- 0.00011*(1.0-thk*thk) / 0.875    ) / fc;
	cf_rt  = gex*cfo/(fc*grt) / rt;
	cf_msq = gex*cfo/(fc*grt) * (-0.25*gm1/fc/fc) - 0.25*gm1*(cf)/fc/fc;
	
% 	return true;
% }
end
function [di, di_hk, di_rt]=dil(hk, rt)
% ////--------------------------------------
% //     dil function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% //---- laminar dissipation function  ( 2 cd/h* )     (from falkner-skan)
	if(hk<4.0) 
		di    = ( 0.00205   *  ((4.0-hk)^5.5) + 0.207 ) / rt;
		di_hk = ( -.00205*5.5*((4.0-hk)^4.5)         ) / rt;
	
    else
		hkb = hk - 4.0;
		den = 1.0 + 0.02*hkb*hkb;
		di    = ( -.0016  *  hkb*hkb  /den   + 0.207              ) / rt;
		di_hk = ( -.0016*2.0*hkb*(1.0/den - 0.02*hkb*hkb/den/den) ) / rt;
    end
	di_rt = -(di)/rt;
	
% 	return true;
% }
end
function [di, di_hk, di_rt]=dilw( hk,  rt)
% ////--------------------------------------
% //     dilw function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% 	//	double msq = 0.0;
% 	double hs, hs_hk, hs_rt, hs_msq;

	[hs, hs_hk, hs_rt, hs_msq]=hsl(hk);
% 	//---- laminar wake dissipation function  ( 2 cd/h* )
	rcd    =  1.10 * (1.0 - 1.0/hk)* (1.0 - 1.0/hk) / hk;
	rcd_hk = -1.10 * (1.0 - 1.0/hk)*2.0/hk/hk/hk- rcd/hk;
	
	di    = 2.0*rcd   /(hs*rt);
	di_hk = 2.0*rcd_hk/(hs*rt) - ((di)/hs)*hs_hk;
	di_rt = -(di)/rt         - ((di)/hs)*hs_rt;
	
% 	return true;
% }
end
function [p]=blkin(p)
% ////--------------------------------------
% //     blkin function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

%//----------------------------------------------------------
%//     calculates turbulence-independent secondary "2" 
%//     variables from the primary "2" variables.
%//----------------------------------------------------------
% 	double p.tr2, p.herat, p.he_u2, p.he_ms, v2_he, p.hk2_h2, p.hk2_m2;
	%//---- set edge mach number ** 2
	p.m2    = p.u2*p.u2*p.hstinv / (p.gm1bl*(1.0 - 0.5*p.u2*p.u2*p.hstinv));
	p.tr2   = 1.0 + 0.5*p.gm1bl*p.m2;
	p.m2_u2 = 2.0*p.m2*p.tr2/p.u2;
	p.m2_ms = p.u2*p.u2*p.tr2 / (p.gm1bl*(1.0 - 0.5*p.u2*p.u2*p.hstinv))* p.hstinv_ms;
	
	%//---- set edge static density (isentropic relation)
	p.r2    = p.rstbl   *(p.tr2^(-1.0/p.gm1bl));
	p.r2_u2 = -p.r2/p.tr2 * 0.5*p.m2_u2;
	p.r2_ms = -p.r2/p.tr2 * 0.5*p.m2_ms+ p.rstbl_ms*(p.tr2^(-1.0/p.gm1bl));
	
	%//---- set shape parameter
	p.h2    =  p.d2/p.t2;
	p.h2_d2 = 1.0/p.t2;
	p.h2_t2 = -p.h2/p.t2;
	
	%//---- set edge static/stagnation enthalpy
	p.herat = 1.0 - 0.5*p.u2*p.u2*p.hstinv;
	p.he_u2 =      -        p.u2*p.hstinv;
	p.he_ms =      - 0.5*p.u2*p.u2*p.hstinv_ms;
	%//---- set molecular viscosity
	p.v2 = sqrt(p.herat*p.herat*p.herat) * (1.0+p.hvrat)/(p.herat+p.hvrat)/p.reybl;
	v2_he = p.v2*(1.5/p.herat - 1.0/(p.herat+p.hvrat));
	
	p.v2_u2 =                        v2_he*p.he_u2;
	p.v2_ms = -p.v2/p.reybl * p.reybl_ms + v2_he*p.he_ms;
	p.v2_re = -p.v2/p.reybl * p.reybl_re;
	
	%//---- set kinematic shape parameter
	[ p.hk2, p.hk2_h2, p.hk2_m2]=hkin(p.h2, p.m2 );
	
	p.hk2_u2 =                p.hk2_m2*p.m2_u2;
	p.hk2_t2 = p.hk2_h2*p.h2_t2;
	p.hk2_d2 = p.hk2_h2*p.h2_d2;
	p.hk2_ms =                p.hk2_m2*p.m2_ms;
	
	%//---- set momentum thickness reynolds number
	p.rt2    = p.r2*p.u2*p.t2/p.v2;
	p.rt2_u2 = p.rt2*(1.0/p.u2 + p.r2_u2/p.r2 - p.v2_u2/p.v2);
	p.rt2_t2 = p.rt2/p.t2;
	p.rt2_ms = p.rt2*(         p.r2_ms/p.r2 - p.v2_ms/p.v2);
	p.rt2_re = p.rt2*(                  - p.v2_re/p.v2);
	
% 	return true;
% }
end
function [hk,hk_h,hk_msq]=hkin(h, msq)
% ////--------------------------------------
% //     hkin function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% 	//---- calculate kinematic shape parameter (assuming air)
% 	//     (from Whitfield )
	hk     =    (h - 0.29*msq)   /(1.0 + 0.113*msq);
	hk_h   =     1.0              /(1.0 + 0.113*msq);
	hk_msq = (-.29 - 0.113*(hk))/(1.0 + 0.113*msq);
	
% 	return true;
% }
end
function [p]=blprv(xsi, ami, cti, thi, dsi,  dswaki,  uei,p)
% ////--------------------------------------
% //     blprv function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% //----------------------------------------------------------
% //     set bl primary "2" variables from parameter list
% //----------------------------------------------------------

	p.x2 = xsi;
	p.ampl2 = ami;
	p.s2  = cti;
	p.t2  = thi;
	p.d2  = dsi - dswaki;
	p.dw2 = dswaki;
	
	p.u2 =uei*(1.0-p.tkbl) / (1.0 - p.tkbl*(uei/p.qinfbl)*(uei/p.qinfbl));
	p.u2_uei = (1.0 + p.tkbl*(2.0*p.u2*uei/p.qinfbl/p.qinfbl - 1.0))...
		/ (1.0 - p.tkbl*(uei/p.qinfbl)*(uei/p.qinfbl));
	p.u2_ms  = (p.u2*(uei/p.qinfbl)*(uei/p.qinfbl)  -  uei)*p.tkbl_ms...
		/ (1.0 - p.tkbl*(uei/p.qinfbl)*(uei/p.qinfbl));
% 	return true;

% }
end
function [p]=xifset(is,p)
% ////--------------------------------------
% //     xifset function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% 	%//-----------------------------------------------------
% 	%//	   sets forced-transition bl coordinate locations.
% 	%//-----------------------------------------------------
% 	double chx, chy, chsq, str;

	if(p.xstrip(is)>=1.0) 
		p.xiforc = p.xssi(p.iblte(is),is);
		return ;
    end
	
	p.chx = p.xte - p.xle;
	p.chy = p.yte - p.yle;
	p.chsq = p.chx*p.chx + p.chy*p.chy;
	
	%//---- calculate chord-based x/c, y/c
	for(i=1: p.n)
		p.w1(i) = ((p.x(i)-p.xle)*p.chx + (p.y(i)-p.yle)*p.chy) / p.chsq;
		p.w2(i) = ((p.y(i)-p.yle)*p.chx - (p.x(i)-p.xle)*p.chy) / p.chsq;
    end
	
	p.w3=SPLIND(p.w1,p.s,p.n,-999.0,-999.0);
	p.w4=SPLIND(p.w2,p.s,p.n,-999.0,-999.0);
	
	if(is==1) 
		
		%//----- set approximate arc length of forced transition point for sinvrt
		p.str = p.sle + (p.s(1)-p.sle)*p.xstrip(is);
		
		%//----- calculate actual arc length
		p.str=sinvrt(p.str,p.xstrip(is),p.w1,p.w3,p.s,p.n);
		
		%//----- set bl coordinate value
		p.xiforc = min((p.sst-p.str), p.xssi(p.iblte(is),is));
	
    else
		%//----- same for bottom side
		
		p.str = p.sle + (p.s(n)-p.sle)*p.xstrip(is);
		p.str=sinvrt(p.str,p.xstrip(is),p.w1,p.w3,p.s,p.n);
		p.xiforc = min((p.str - p.sst) , p.xssi(p.iblte(is),is));
		
    end
	
	if(p.xiforc < 0.0) 
		%//TRACE(" ***  stagnation point is past trip on side %d\n", is);
% 		CString str;
		disp(' ***  stagnation point is past trip on side ');
% 		if(m_bTrace)pXFile->WriteString(str);

		p.xiforc = p.xssi(p.iblte(is),is);
    end
	
% 	return true;
% }
end
function [p]=blsys(p)
% ////--------------------------------------
% //     blsys function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

	%//------------------------------------------------------------------
	%//
	%//     sets up the bl newton system governing the current interval:
	%//
	%//     |       ||da1|     |       ||da2|       |     |
	%//     |  vs1  ||dt1|  +  |  vs2  ||dt2|   =   |vsrez|
	%//     |       ||dd1|     |       ||dd2|       |     |
	%//              |du1|              |du2|
	%//              |dx1|              |dx2|
	%//
	%//        3x5    5x1         3x5    5x1          3x1
	%//
	%//     the system as shown corresponds to a laminar station
	%//     if tran, then  ds2  replaces  da2
	%//     if turb, then  ds1, ds2  replace  da1, da2
	%//
	%//------------------------------------------------------------------
	%//      implicit real(m)
% 	double res_u1, res_u2, res_ms;
	
	%//---- calculate secondary bl variables and their sensitivities
    p.bule=1;
	if(p.wake) 
        [p]=blvar(3,p);
% 		blvar(3);
% 		blmid(3);
        [p]=blmid(3,p);
    else 
		if(p.turb || p.tran) 
% 			blvar(2);
            [p]=blvar(2,p);
            % 			blmid(2);
     
            [p]=blmid(2,p);		
        else
% 			blvar(1);
            [p]=blvar(1,p);
            % 			blmid(1);
     
            [p]=blmid(1,p);	
        end
    end

	%//---- for the similarity station, "1" and "2" variables are the same
	if(p.simi) 
		%//		for(int icom=1;icom<= ncom;icom++) com1(icom) = com2(icom);
		[p]=stepbl(p);
		
    end

	%//---- set up appropriate finite difference system for current interval
	if(p.tran)

		p=trdif(p);

    elseif(p.simi)
% 		bldif(0);
       [p]=bldif(0,p);
    elseif(p.turb==false)
% 		bldif(1);
[p]=bldif(1,p);
    elseif(p.wake)
% 		bldif(3);
  [p]=bldif(3,p);
    elseif(p.turb)
% 		bldif(2);
[p]=bldif(2,p);
    end

	if(p.simi) 
		%//----- at similarity station, "1" variables are really "2" variables
		for ( k=1: 4)
			for( l=1: 5)
				p.vs2(k,l) = p.vs1(k,l) + p.vs2(k,l);
				p.vs1(k,l) = 0.0;
            end
        end
    end
	
	%//---- change system over into incompressible uei and mach
	for(k=1: 4)
		
        
		%//------ residual derivatives wrt compressible uec
		p.res_u1 = p.vs1(k,4);
		p.res_u2 = p.vs2(k,4);
		p.res_ms = p.vsm(k);
		
		%//------ combine with derivatives of compressible  u1,u2 = uec(uei m)
		p.vs1(k,4) = p.res_u1*p.u1_uei;
		p.vs2(k,4) =                p.res_u2*p.u2_uei;
		p.vsm(k)   = p.res_u1*p.u1_ms + p.res_u2*p.u2_ms  + p.res_ms;
    end
% 	return true;
% 
end
function [p]=stepbl(p)
% ////--------------------------------------
% //     stepbl function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% 	//arcds : can't think of a more elegant way to do this, and too lazy to search
p.x1=p.x2;
p.u1=p.u2;
p.t1=p.t2;
p.d1=p.d2;
p.s1=p.s2;
p.ampl1  =p.ampl2;
p.u1_uei =p.u2_uei;
p.u1_ms  =p.u2_ms;
p.dw1=p.dw2;
p.h1 =p.h2;
p.h1_t1  =p.h2_t2;
p.h1_d1  =p.h2_d2;
p.m1=p.m2;
p.m1_u1  =p.m2_u2;
p.m1_ms  =p.m2_ms;
p.r1=p.r2;
p.r1_u1  =p.r2_u2;
p.r1_ms  =p.r2_ms;
p.v1=p.v2;
p.v1_u1  =p.v2_u2;
p.v1_ms  =p.v2_ms;
p.v1_re  =p.v2_re;
p.hk1=p.hk2;
p.hk1_u1 =p.hk2_u2;
p.hk1_t1 =p.hk2_t2;
p.hk1_d1 =p.hk2_d2;
p.hk1_ms =p.hk2_ms;
p.hs1=p.hs2;
p.hs1_u1 =p.hs2_u2;
p.hs1_t1 =p.hs2_t2;
p.hs1_d1 =p.hs2_d2;
p.hs1_ms =p.hs2_ms;
p.hs1_re =p.hs2_re;
p.hc1=p.hc2;
p.hc1_u1 =p.hc2_u2;
p.hc1_t1 =p.hc2_t2;
p.hc1_d1 =p.hc2_d2;
p.hc1_ms =p.hc2_ms;
p.rt1=p.rt2;
p.rt1_u1 =p.rt2_u2;
p.rt1_t1 =p.rt2_t2;
p.rt1_ms =p.rt2_ms;
p.rt1_re =p.rt2_re;
p.cf1=p.cf2;
p.cf1_u1 =p.cf2_u2;
p.cf1_t1 =p.cf2_t2;
p.cf1_d1 =p.cf2_d2;
p.cf1_ms =p.cf2_ms;
p.cf1_re =p.cf2_re;
p.di1=p.di2;
p.di1_u1 =p.di2_u2;
p.di1_t1 =p.di2_t2;
p.di1_d1 =p.di2_d2;
p.di1_s1 =p.di2_s2;
p.di1_ms =p.di2_ms;
p.di1_re =p.di2_re;
p.us1=p.us2;
p.us1_u1 =p.us2_u2;
p.us1_t1 =p.us2_t2;
p.us1_d1 =p.us2_d2;
p.us1_ms =p.us2_ms;
p.us1_re =p.us2_re;
p.cq1=p.cq2;
p.cq1_u1 =p.cq2_u2;
p.cq1_t1 =p.cq2_t2;
p.cq1_d1 =p.cq2_d2;
p.cq1_ms =p.cq2_ms;
p.cq1_re =p.cq2_re;
p.de1=p.de2;
p.de1_u1 =p.de2_u2;
p.de1_t1 =p.de2_t2;
p.de1_d1 =p.de2_d2;
p.de1_ms =p.de2_ms;
% 	return true;
end
function p=ueset(p)
% ////--------------------------------------
% //     ueset function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% 	//---------------------------------------------------------
% 	//     sets ue from inviscid ue plus all source influence
% 	//---------------------------------------------------------
% 	int i, is,ibl, j, js, jbl;
% 	double dui, ue_m;
	for (is=1:2)
		for(ibl=2: p.nbl(is))
			i = p.ipan(ibl,is);
			
			p.dui = 0.0;
			for (js=1: 2)
				for(jbl=2: p.nbl(js))
					j  = p.ipan(jbl,js);
					p.ue_m = -p.vti(ibl,is)*p.vti(jbl,js)*p.dij(i,j);
					p.dui = p.dui + p.ue_m*p.mass(jbl,js);
                end
            end
			
			p.uedg(ibl,is) = p.uinv(ibl,is) + p.dui;
			
        end
    end
% 	return true;
% }
end
function p=trdif(p)
% ////--------------------------------------
% //     trdif function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------


%//-----------------------------------------------
%//     sets up the newton system governing the
%//     transition interval.  equations governing
%//     the  laminar  part  p.x1 < xi < p.xt  and
%//     the turbulent part  p.xt < xi < p.x2
%//     are simply summed.
%//-----------------------------------------------
%     double bl1(5,6), bl2(5,6), blrez(5), blm(5), blr(5), blx(5),
% 		   p.bt1(5,6), p.bt2(5,6), p.btrez(5), btm(5), btr(5), btx(5);
% 	double p.wf2, p.wf2_xt, p.wf2_a1,p.wf2_x1, p.wf2_x2, p.wf2_t1, p.wf2_t2;
% 	double p.wf2_d1, p.wf2_d2, p.wf2_u1, p.wf2_u2, p.wf2_ms, p.wf2_re, p.wf2_xf;
% 	double p.wf1, p.wf1_a1, p.wf1_x1, p.wf1_x2, p.wf1_t1, p.wf1_t2, p.wf1_d1, p.wf1_d2;
% 	double p.wf1_u1, p.wf1_u2, p.wf1_ms, p.wf1_re, p.wf1_xf;
% 	double p.tt, p.tt_a1, p.tt_x1, p.tt_x2, p.tt_t1, p.tt_t2, p.tt_d1, p.tt_d2, p.tt_u1, p.tt_u2;
% 	double p.tt_ms, p.tt_re, p.tt_xf, p.dt, p.dt_a1, p.dt_x1, p.dt_x2, p.dt_t1, p.dt_t2;
% 	double p.dt_d1, p.dt_d2, p.dt_u1, p.dt_u2, p.dt_ms, p.dt_re, p.dt_xf;
% 	double ut, p.ut_a1, p.ut_x1, p.ut_x2, p.ut_t1, p.ut_t2, p.ut_d1, p.ut_d2, p.ut_u1, p.ut_u2;
% 	double p.ut_ms, p.ut_re, p.ut_xf;
% 	double st, p.st_p.tt, p.st_p.dt, p.st_ut, p.st_ms, p.st_re, p.st_a1, p.st_x1, p.st_x2, p.st_t1, p.st_t2;
% 	double p.st_d1, p.st_d2, p.st_u1, p.st_u2, p.st_xf;
% 	double p.ctr, p.ctr_hk2;
% 	 k;
%//	double c1sav(74), c2sav(74);
	
	%//---- save variables and sensitivities for future restoration
%//	for ( icom=1; icom<= ncom;icom++){
%//		c1sav(icom) = com1(icom);
%//		c2sav(icom) = com2(icom);
%//	}
     blsav(1).d=0;
 	[blsav,p]=saveblData(1,p,blsav);

 	[blsav,p]=saveblData(2,p,blsav);
    
	
	%//---- weighting factors for linear interpolation to transition point
	p.wf2    = (p.xt-p.x1)/(p.x2-p.x1);
	p.wf2_xt = 1.0/(p.x2-p.x1);
	
	p.wf2_a1 = p.wf2_xt*p.xt_a1;
	p.wf2_x1 = p.wf2_xt*p.xt_x1 + (p.wf2-1.0)/(p.x2-p.x1);
	p.wf2_x2 = p.wf2_xt*p.xt_x2 -  p.wf2      /(p.x2-p.x1);
	p.wf2_t1 = p.wf2_xt*p.xt_t1;
	p.wf2_t2 = p.wf2_xt*p.xt_t2;
	p.wf2_d1 = p.wf2_xt*p.xt_d1;
	p.wf2_d2 = p.wf2_xt*p.xt_d2;
	p.wf2_u1 = p.wf2_xt*p.xt_u1;
	p.wf2_u2 = p.wf2_xt*p.xt_u2;
	p.wf2_ms = p.wf2_xt*p.xt_ms;
	p.wf2_re = p.wf2_xt*p.xt_re;
	p.wf2_xf = p.wf2_xt*p.xt_xf;
	
	p.wf1    = 1.0 - p.wf2;
	p.wf1_a1 = -p.wf2_a1;
	p.wf1_x1 = -p.wf2_x1;
	p.wf1_x2 = -p.wf2_x2;
	p.wf1_t1 = -p.wf2_t1;
	p.wf1_t2 = -p.wf2_t2;
	p.wf1_d1 = -p.wf2_d1;
	p.wf1_d2 = -p.wf2_d2;
	p.wf1_u1 = -p.wf2_u1;
	p.wf1_u2 = -p.wf2_u2;
	p.wf1_ms = -p.wf2_ms;
	p.wf1_re = -p.wf2_re;
	p.wf1_xf = -p.wf2_xf;
	
	%//**** first,  do laminar part between p.x1 and p.xt
	
	%//-----interpolate primary variables to transition point
	p.tt    = p.t1*p.wf1    + p.t2*p.wf2;
	p.tt_a1 = p.t1*p.wf1_a1 + p.t2*p.wf2_a1;
	p.tt_x1 = p.t1*p.wf1_x1 + p.t2*p.wf2_x1;
	p.tt_x2 = p.t1*p.wf1_x2 + p.t2*p.wf2_x2;
	p.tt_t1 = p.t1*p.wf1_t1 + p.t2*p.wf2_t1 + p.wf1;
	p.tt_t2 = p.t1*p.wf1_t2 + p.t2*p.wf2_t2 + p.wf2;
	p.tt_d1 = p.t1*p.wf1_d1 + p.t2*p.wf2_d1;
	p.tt_d2 = p.t1*p.wf1_d2 + p.t2*p.wf2_d2;
	p.tt_u1 = p.t1*p.wf1_u1 + p.t2*p.wf2_u1;
	p.tt_u2 = p.t1*p.wf1_u2 + p.t2*p.wf2_u2;
	p.tt_ms = p.t1*p.wf1_ms + p.t2*p.wf2_ms;
	p.tt_re = p.t1*p.wf1_re + p.t2*p.wf2_re;
	p.tt_xf = p.t1*p.wf1_xf + p.t2*p.wf2_xf;
	
	p.dt    = p.d1*p.wf1    + p.d2*p.wf2;
	p.dt_a1 = p.d1*p.wf1_a1 + p.d2*p.wf2_a1;
	p.dt_x1 = p.d1*p.wf1_x1 + p.d2*p.wf2_x1;
	p.dt_x2 = p.d1*p.wf1_x2 + p.d2*p.wf2_x2;
	p.dt_t1 = p.d1*p.wf1_t1 + p.d2*p.wf2_t1;
	p.dt_t2 = p.d1*p.wf1_t2 + p.d2*p.wf2_t2;
	p.dt_d1 = p.d1*p.wf1_d1 + p.d2*p.wf2_d1 + p.wf1;
	p.dt_d2 = p.d1*p.wf1_d2 + p.d2*p.wf2_d2 + p.wf2;
	p.dt_u1 = p.d1*p.wf1_u1 + p.d2*p.wf2_u1;
	p.dt_u2 = p.d1*p.wf1_u2 + p.d2*p.wf2_u2;
	p.dt_ms = p.d1*p.wf1_ms + p.d2*p.wf2_ms;
	p.dt_re = p.d1*p.wf1_re + p.d2*p.wf2_re;
	p.dt_xf = p.d1*p.wf1_xf + p.d2*p.wf2_xf;
	
	p.ut    = p.u1*p.wf1    + p.u2*p.wf2;
	p.ut_a1 = p.u1*p.wf1_a1 + p.u2*p.wf2_a1;
	p.ut_x1 = p.u1*p.wf1_x1 + p.u2*p.wf2_x1;
	p.ut_x2 = p.u1*p.wf1_x2 + p.u2*p.wf2_x2;
	p.ut_t1 = p.u1*p.wf1_t1 + p.u2*p.wf2_t1;
	p.ut_t2 = p.u1*p.wf1_t2 + p.u2*p.wf2_t2;
	p.ut_d1 = p.u1*p.wf1_d1 + p.u2*p.wf2_d1;
	p.ut_d2 = p.u1*p.wf1_d2 + p.u2*p.wf2_d2;
	p.ut_u1 = p.u1*p.wf1_u1 + p.u2*p.wf2_u1 + p.wf1;
	p.ut_u2 = p.u1*p.wf1_u2 + p.u2*p.wf2_u2 + p.wf2;
	p.ut_ms = p.u1*p.wf1_ms + p.u2*p.wf2_ms;
	p.ut_re = p.u1*p.wf1_re + p.u2*p.wf2_re;
	p.ut_xf = p.u1*p.wf1_xf + p.u2*p.wf2_xf;
	
	%//---- set primary "t" variables at p.xt  (really placed into "2" variables)
	p.x2 = p.xt;
	p.t2 = p.tt;
	p.d2 = p.dt;
	p.u2 = p.ut;

	p.ampl2 = p.amcrit;
	p.s2 = 0.0;
	
	%//---- calculate laminar secondary "t" variables
	p=blkin(p);
	p=blvar(1,p);
	
	%//---- calculate p.x1-p.xt midpoint cfm value
	p=blmid(1,p);
	
	%//=    at this point, all "2" variables are really "t" variables at p.xt
	
	
	%//---- set up newton system for dam, p.dth, dds, due, dxi  at  p.x1 and p.xt
	p=bldif(1,p);
	
	%//---- the current newton system is in terms of "1" and "t" variables,
	%//-    so calculate its equivalent in terms of "1" and "2" variables.
	%//-    in other words, convert residual sensitivities wrt "t" variables
	%//-    into sensitivities wrt "1" and "2" variables.  the amplification
	%//-    equation is unnecessary here, so the k=1 row is left empty.
	for (k=2: 3)
		p.blrez(k) = p.vsrez(k);
		p.blm(k)   = p.vsm(k)+ p.vs2(k,2)*p.tt_ms+ p.vs2(k,3)*p.dt_ms+ p.vs2(k,4)*p.ut_ms+ p.vs2(k,5)*p.xt_ms;
		p.blr(k)   = p.vsr(k)+ p.vs2(k,2)*p.tt_re+ p.vs2(k,3)*p.dt_re+ p.vs2(k,4)*p.ut_re+ p.vs2(k,5)*p.xt_re;
		p.blx(k)   = p.vsx(k)+ p.vs2(k,2)*p.tt_xf+ p.vs2(k,3)*p.dt_xf+ p.vs2(k,4)*p.ut_xf+ p.vs2(k,5)*p.xt_xf;
		
		p.bl1(k,1) = p.vs1(k,1)+ p.vs2(k,2)*p.tt_a1+ p.vs2(k,3)*p.dt_a1+ p.vs2(k,4)*p.ut_a1+ p.vs2(k,5)*p.xt_a1;
		p.bl1(k,2) = p.vs1(k,2)+ p.vs2(k,2)*p.tt_t1+ p.vs2(k,3)*p.dt_t1+ p.vs2(k,4)*p.ut_t1+ p.vs2(k,5)*p.xt_t1;
		p.bl1(k,3) = p.vs1(k,3)+ p.vs2(k,2)*p.tt_d1+ p.vs2(k,3)*p.dt_d1+ p.vs2(k,4)*p.ut_d1+ p.vs2(k,5)*p.xt_d1;
		p.bl1(k,4) = p.vs1(k,4)+ p.vs2(k,2)*p.tt_u1+ p.vs2(k,3)*p.dt_u1+ p.vs2(k,4)*p.ut_u1+ p.vs2(k,5)*p.xt_u1;
		p.bl1(k,5) = p.vs1(k,5)+ p.vs2(k,2)*p.tt_x1+ p.vs2(k,3)*p.dt_x1+ p.vs2(k,4)*p.ut_x1+ p.vs2(k,5)*p.xt_x1;
		
		p.bl2(k,1) = 0.0;
		p.bl2(k,2) = p.vs2(k,2)*p.tt_t2+ p.vs2(k,3)*p.dt_t2+ p.vs2(k,4)*p.ut_t2+ p.vs2(k,5)*p.xt_t2;
		p.bl2(k,3) = p.vs2(k,2)*p.tt_d2+ p.vs2(k,3)*p.dt_d2+ p.vs2(k,4)*p.ut_d2+ p.vs2(k,5)*p.xt_d2;
		p.bl2(k,4) = p.vs2(k,2)*p.tt_u2+ p.vs2(k,3)*p.dt_u2+ p.vs2(k,4)*p.ut_u2+ p.vs2(k,5)*p.xt_u2;
		p.bl2(k,5) = p.vs2(k,2)*p.tt_x2+ p.vs2(k,3)*p.dt_x2+ p.vs2(k,4)*p.ut_x2+ p.vs2(k,5)*p.xt_x2;
		
    end
	
	%//**** second, set up turbulent part between p.xt and p.x2  ****
	
	%//---- calculate equilibrium shear coefficient cqt at transition point
	p=blvar(2,p);
	
	%//---- set initial shear coefficient value st at transition point
	%//-    ( note that cq2, cq2_t2, etc. are really "cqt", "cqt_p.tt", etc.)
	
	p.ctr     = 1.8*exp(-3.3/(p.hk2-1.0));
	p.ctr_hk2 = p.ctr * 3.3/(p.hk2-1.0)/(p.hk2-1.0);
	
	p.st    = p.ctr*p.cq2;
	p.st_tt = p.ctr*p.cq2_t2 + p.cq2*p.ctr_hk2*p.hk2_t2;
	p.st_p.dt = p.ctr*p.cq2_d2 + p.cq2*p.ctr_hk2*p.hk2_d2;
	p.st_ut = p.ctr*p.cq2_u2 + p.cq2*p.ctr_hk2*p.hk2_u2;
	p.st_ms = p.ctr*p.cq2_ms + p.cq2*p.ctr_hk2*p.hk2_ms;
	p.st_re = p.ctr*p.cq2_re;
	
	%//---- calculate st sensitivities wrt the actual "1" and "2" variables
	p.st_a1 = p.st_tt*p.tt_a1 + p.st_p.dt*p.dt_a1 + p.st_ut*p.ut_a1;
	p.st_x1 = p.st_tt*p.tt_x1 + p.st_p.dt*p.dt_x1 + p.st_ut*p.ut_x1;
	p.st_x2 = p.st_tt*p.tt_x2 + p.st_p.dt*p.dt_x2 + p.st_ut*p.ut_x2;
	p.st_t1 = p.st_tt*p.tt_t1 + p.st_p.dt*p.dt_t1 + p.st_ut*p.ut_t1;
	p.st_t2 = p.st_tt*p.tt_t2 + p.st_p.dt*p.dt_t2 + p.st_ut*p.ut_t2;
	p.st_d1 = p.st_tt*p.tt_d1 + p.st_p.dt*p.dt_d1 + p.st_ut*p.ut_d1;
	p.st_d2 = p.st_tt*p.tt_d2 + p.st_p.dt*p.dt_d2 + p.st_ut*p.ut_d2;
	p.st_u1 = p.st_tt*p.tt_u1 + p.st_p.dt*p.dt_u1 + p.st_ut*p.ut_u1;
	p.st_u2 = p.st_tt*p.tt_u2 + p.st_p.dt*p.dt_u2 + p.st_ut*p.ut_u2;
	p.st_ms = p.st_tt*p.tt_ms + p.st_p.dt*p.dt_ms + p.st_ut*p.ut_ms + p.st_ms;
	p.st_re = p.st_tt*p.tt_re + p.st_p.dt*p.dt_re + p.st_ut*p.ut_re + p.st_re;
	p.st_xf = p.st_tt*p.tt_xf + p.st_p.dt*p.dt_xf + p.st_ut*p.ut_xf;
	
	p.ampl2 = 0.0;
	p.s2 = p.st;
	
	%//---- recalculate turbulent secondary "t" variables using proper cti
	p=blvar(2,p);
	
	%//---- set "1" variables to "t" variables and reset "2" variables
	%//-    to their saved turbulent values
%//	for (icom=1; icom<= ncom; icom++){
%//		com1(icom) = com2(icom);
%//		com2(icom) = c2sav(icom);
%//	}
    
	p=stepbl(p);
 	[p,blsav]=restoreblData(2,blsav,p);

	
	%//---- calculate p.xt-p.x2 midpoint cfm value
	p=blmid(2,p);
	
	%//---- set up newton system for dct, p.dth, dds, due, dxi  at  p.xt and p.x2
	p=bldif(2,p);
	
	%//---- convert sensitivities wrt "t" variables into sensitivities
	%//-    wrt "1" and "2" variables as done before for the laminar part
	for (k=1:3)
		p.btrez(k) = p.vsrez(k);
		p.btm(k)   = p.vsm(k) + p.vs1(k,1)*p.st_ms+ p.vs1(k,2)*p.tt_ms+ p.vs1(k,3)*p.dt_ms+ p.vs1(k,4)*p.ut_ms+ p.vs1(k,5)*p.xt_ms;
		p.btr(k)   = p.vsr(k) + p.vs1(k,1)*p.st_re+ p.vs1(k,2)*p.tt_re+ p.vs1(k,3)*p.dt_re+ p.vs1(k,4)*p.ut_re+ p.vs1(k,5)*p.xt_re;
		p.btx(k)   = p.vsx(k) + p.vs1(k,1)*p.st_xf+ p.vs1(k,2)*p.tt_xf+ p.vs1(k,3)*p.dt_xf+ p.vs1(k,4)*p.ut_xf+ p.vs1(k,5)*p.xt_xf;
		
		p.bt1(k,1) = p.vs1(k,1)*p.st_a1+ p.vs1(k,2)*p.tt_a1+ p.vs1(k,3)*p.dt_a1+ p.vs1(k,4)*p.ut_a1+ p.vs1(k,5)*p.xt_a1;
		p.bt1(k,2) = p.vs1(k,1)*p.st_t1+ p.vs1(k,2)*p.tt_t1+ p.vs1(k,3)*p.dt_t1+ p.vs1(k,4)*p.ut_t1+ p.vs1(k,5)*p.xt_t1;
		p.bt1(k,3) = p.vs1(k,1)*p.st_d1+ p.vs1(k,2)*p.tt_d1+ p.vs1(k,3)*p.dt_d1+ p.vs1(k,4)*p.ut_d1+ p.vs1(k,5)*p.xt_d1;
		p.bt1(k,4) = p.vs1(k,1)*p.st_u1+ p.vs1(k,2)*p.tt_u1+ p.vs1(k,3)*p.dt_u1+ p.vs1(k,4)*p.ut_u1+ p.vs1(k,5)*p.xt_u1;
		p.bt1(k,5) = p.vs1(k,1)*p.st_x1+ p.vs1(k,2)*p.tt_x1+ p.vs1(k,3)*p.dt_x1+ p.vs1(k,4)*p.ut_x1+ p.vs1(k,5)*p.xt_x1;
		
		p.bt2(k,1) = p.vs2(k,1);
		p.bt2(k,2) = p.vs2(k,2)+ p.vs1(k,1)*p.st_t2+ p.vs1(k,2)*p.tt_t2+ p.vs1(k,3)*p.dt_t2+ p.vs1(k,4)*p.ut_t2+ p.vs1(k,5)*p.xt_t2;
		p.bt2(k,3) = p.vs2(k,3)+ p.vs1(k,1)*p.st_d2+ p.vs1(k,2)*p.tt_d2+ p.vs1(k,3)*p.dt_d2+ p.vs1(k,4)*p.ut_d2+ p.vs1(k,5)*p.xt_d2;
		p.bt2(k,4) = p.vs2(k,4)+ p.vs1(k,1)*p.st_u2+ p.vs1(k,2)*p.tt_u2+ p.vs1(k,3)*p.dt_u2+ p.vs1(k,4)*p.ut_u2+ p.vs1(k,5)*p.xt_u2;
		p.bt2(k,5) = p.vs2(k,5)+ p.vs1(k,1)*p.st_x2+ p.vs1(k,2)*p.tt_x2+ p.vs1(k,3)*p.dt_x2+ p.vs1(k,4)*p.ut_x2+ p.vs1(k,5)*p.xt_x2;
		
    end
	
	%//---- add up laminar and turbulent parts to get final system
	%//-    in terms of honest-to-god "1" and "2" variables.
	p.vsrez(1) =            p.btrez(1);
	p.vsrez(2) = p.blrez(2) + p.btrez(2);
	p.vsrez(3) = p.blrez(3) + p.btrez(3);
	p.vsm(1)   =            p.btm(1);
	p.vsm(2)   = p.blm(2)   + p.btm(2);
	p.vsm(3)   = p.blm(3)   + p.btm(3);
	p.vsr(1)   =            p.btr(1);
	p.vsr(2)   = p.blr(2)   + p.btr(2);
	p.vsr(3)   = p.blr(3)   + p.btr(3);
	p.vsx(1)   =            p.btx(1);
	p.vsx(2)   = p.blx(2)   + p.btx(2);
	p.vsx(3)   = p.blx(3)   + p.btx(3);
	for ( l=1:5)
		p.vs1(1,l) =            p.bt1(1,l);
		p.vs2(1,l) =            p.bt2(1,l);
		p.vs1(2,l) = p.bl1(2,l) + p.bt1(2,l);
		p.vs2(2,l) = p.bl2(2,l) + p.bt2(2,l);
		p.vs1(3,l) = p.bl1(3,l) + p.bt1(3,l);
		p.vs2(3,l) = p.bl2(3,l) + p.bt2(3,l);
    end
	
	%//---- to be sanitary, restore "1" quantities which got clobbered
	%//-    in all of the numerical gymnastics above.  the "2" variables
	%//-    were already restored for the p.xt-p.x2 differencing part.
%//	for (icom=1; icom<=ncom;icom++){
%//		com1(icom) = c1sav(icom);
%//	}
 	[p,blsav]=restoreblData(1,blsav,p);
	
% 	return true;		
% }
end
function p=mrchdu(p)
% ////--------------------------------------
% //     mrchdu function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% {
% 	%//----------------------------------------------------
% 	%//     marches the bls and p.wake in mixed mode using
% 	%//     the current ue and hk.  the calculated ue
% 	%//     and hk lie along a line quasi-normal to the
% 	%//     natural ue-hk charap.cteristic line of the
% 	%//     current bl so that the goldstein or levy-lees
% 	%//     singularity is never encountered.  continuous
% 	%//     checking of transition onset is performed.
% 	%//----------------------------------------------------
% 	CString str;
% 
% 	vtmp(5,6), vztmp(5);
% 	
 	deps = 0.000005;
% 	int ibl, ibm, p.itrold, iw, itbl;%//icom
% 
% 	senswt, thi, p.uei, dsi, p.cti, p.dswaki, ratlen;
% 	sens, p.sennew, dummy, p.msq, p.thm, p.dsm, p.uem;
% 	p.xsi, p.hklim, p.dsw;
	p.ami = 0.0;%//added arcds
	p.dte = 0.0;
	p.cte = 0.0;
	p.tte = 0.0;
	p.ueref = 0.0;
	p.hkref = 0.0;
	p.dmax = 0.0;
	%//---- constant controlling how far hk is allowed to deviate
	%//    from the specified value.
	senswt = 1000.0;
	sens = 0.0;
	p.sennew = 0.0;
	
	for (is=1: 2)%//2000
		
		%//---- set forced transition arc length position
		p=xifset(is,p);
		
		%//---- set leading edge pressure gradient parameter  x/u du/dx
		ibl = 2;
		p.xsi = p.xssi(ibl,is);
		p.uei = p.uedg(ibl,is);
		
		
		%//---- old transition station
		p.itrold = p.itran(is);
		
		p.tran = false;
		p.turb = false;
		p.itran(is) = p.iblte(is);
		%//---- march downstream
		for(ibl=2: p.nbl(is))%//1000
			ibm = ibl-1;

			p.simi = ibl==2;
			p.wake = ibl>p.iblte(is);
			
			%//------ initialize current station to existing variables
			p.xsi = p.xssi(ibl,is);
			p.uei = p.uedg(ibl,is);
			thi = p.thet(ibl,is);
			dsi = p.dstr(ibl,is);
			
			%//------ fixed bug   md 7 june 99
			if(ibl<p.itrold) 
				p.ami = p.ctau(ibl,is);%// p.ami must be initialized
				p.cti = 0.03;
			
            else
				p.cti = p.ctau(ibl,is);
				if(p.cti<=0.0) 
                    p.cti = 0.03;
                end
            end
			
			if(p.wake) 
				iw = ibl - p.iblte(is);
				p.dswaki = p.wgap(iw);
			
            else
                p.dswaki = 0.0;
            end
			
			if(ibl<=p.iblte(is)) 
                dsi = max(dsi-p.dswaki,1.02000*thi) + p.dswaki;
            end
			if(ibl>p.iblte(is))
                dsi = max(dsi-p.dswaki,1.00005*thi) + p.dswaki;
            end

			%//------ newton iteration loop for current station
						
			for (itbl=1: 25)%//100
				
				%//-------- assemble 10x3 linearized system for dctau, dth, dds, due, dxi
				%//         at the previous "1" station and the current "2" station
				%//         (the "1" station coefficients will be ignored)
				
				p=blprv(p.xsi,p.ami,p.cti,thi,dsi,p.dswaki,p.uei,p);
				p=blkin(p);
				
				%//-------- check for transition and set appropriate flags and things
				if((p.simi==false) && (p.turb==false)) 
					p=trchek(p);
					p.ami = p.ampl2;
					if( p.tran) 
                        p.itran(is) = ibl;
                    end
					if(p.tran==false)
                        p.itran(is) = ibl+2;
                    end
                end
				if(ibl==p.iblte(is)+1) 
					p.tte = p.thet(p.iblte(1),1) + p.thet(p.iblte(2),2);
					p.dte = p.dstr(p.iblte(1),1) + p.dstr(p.iblte(2),2) + p.ante;
					p.cte = ( p.ctau(p.iblte(1),1)*p.thet(p.iblte(1),1)...
						+ p.ctau(p.iblte(2),2)*p.thet(p.iblte(2),2) ) / p.tte;
					p=tesys(p.cte,p.tte,p.dte,p);
				
                else
					p=blsys(p);
                end
				
				%//-------- set stuff at first iteration...
				if(itbl==1) 
					
					%//--------- set "baseline" ue and hk for forming  ue(hk)  relation
					p.ueref = p.u2;
					p.hkref = p.hk2;
					
					%//--------- if current point ibl was turbulent and is now laminar, then...
					if(ibl<p.itran(is) && ibl>=p.itrold ) 
						%//---------- extrapolate baseline hk
						p.uem = p.uedg(ibl-1,is);
						p.dsm = p.dstr(ibl-1,is);
						p.thm = p.thet(ibl-1,is);
						p.msq = p.uem*p.uem*p.hstinv / (p.gm1bl*(1.0 - 0.5*p.uem*p.uem*p.hstinv));
						[p.hkref, dummy, dummy]=hkin( p.dsm/p.thm, p.msq );
                    end
					
					%//--------- if current point ibl was laminar, then...
					if(ibl<p.itrold) 
						%//---------- reinitialize or extrapolate p.ctau if it's now turbulent
						if(p.tran)
                            p.ctau(ibl,is) = 0.03;
                        end
						if(p.turb) 
                            p.ctau(ibl,is) = p.ctau(ibl-1,is);
                        end
						if(p.tran || p.turb) 
							p.cti = p.ctau(ibl,is);
							p.s2 = p.cti;
                        end
                    end
                end
				
				if(p.simi || ibl==p.iblte(is)+1) 
					%//--------- for similarity station or first p.wake point, prescribe ue
					p.vs2(4,1) = 0.0;
					p.vs2(4,2) = 0.0;
					p.vs2(4,3) = 0.0;
					p.vs2(4,4) = p.u2_uei;
					p.vsrez(4) = p.ueref - p.u2;
				
                else
					%//******** calculate ue-hk charap.cteristic slope
					for (k=1:4)
						vztmp(k) = p.vsrez(k);
						for (l=1:5)
                            vtmp(k,l) = p.vs2(k,l);
                        end
                    end
					
					
					%//--------- set unit dhk
					vtmp(4,1) = 0.0;
					vtmp(4,2) = p.hk2_t2;
					vtmp(4,3) = p.hk2_d2;
					vtmp(4,4) = p.hk2_u2*p.u2_uei;
					vztmp(4)   = 1.0;
					
					%//--------- calculate due response
					vztmp=Gauss(4,vtmp,vztmp);
					
					%//--------- set  senswt * (normalized due/dhk)
					p.sennew = senswt * vztmp(4) * p.hkref/p.ueref;
					if(itbl<=5)
                        sens = p.sennew;
                    elseif(itbl<=15)
                        sens = 0.5*(sens + p.sennew);
                    end
					
					
					%//--------- set prescribed ue-hk combination
					p.vs2(4,1) = 0.0;
					p.vs2(4,2) =  p.hk2_t2 * p.hkref;
					p.vs2(4,3) =  p.hk2_d2 * p.hkref;
					p.vs2(4,4) =( p.hk2_u2 * p.hkref  +  sens/p.ueref )*p.u2_uei;
					p.vsrez(4)  = -(p.hkref*p.hkref)*(p.hk2 / p.hkref - 1.0)...
								- sens*(p.u2  / p.ueref - 1.0);
					
                end
				
				%//-------- solve newton system for current "2" station
				p.vsrez=Gauss(4,p.vs2,p.vsrez);
				
				%//-------- determine max changes and underrelax if necessary
				p.dmax = max(abs(p.vsrez(2)/thi), abs(p.vsrez(3)/dsi)  );
				if(ibl>=p.itran(is)) 
                    p.dmax = max(p.dmax,abs(p.vsrez(1)/(10.0*p.cti)));
                end
				
				rlx = 1.0;
				if(p.dmax>0.3) 
                    rlx = 0.3/p.dmax;
                end
				
				%//-------- update as usual
				if(ibl<p.itran(is)) 
                    p.ami = p.ami + rlx*p.vsrez(1);
                end
				if(ibl>=p.itran(is)) 
                    p.cti = p.cti + rlx*p.vsrez(1);
                end
				thi = thi + rlx*p.vsrez(2);
				dsi = dsi + rlx*p.vsrez(3);
				p.uei = p.uei + rlx*p.vsrez(4);
				 
				%//-------- eliminate absurd transients
				if(ibl>=p.itran(is)) 
					p.cti = min(p.cti , 0.30);
					p.cti = max(p.cti , 0.0000001);
                end
				
				if(ibl<=p.iblte(is))
                    p.hklim = 1.02;
                else
                    p.hklim = 1.00005;
                end
				
				p.msq = p.uei*p.uei*p.hstinv / (p.gm1bl*(1.0 - 0.5*p.uei*p.uei*p.hstinv));
				p.dsw = dsi - p.dswaki;
				p.dsw=dslim(p.dsw,thi,p.msq,p.hklim);
				dsi = p.dsw + p.dswaki;
				
				if(p.dmax<=deps) 
                    
                    
                    
                    
                    
                    			sens = p.sennew;

			%//------ store primary variables
			if(ibl<p.itran(is)) 
                p.ctau(ibl,is) = p.ami;
            else
                p.ctau(ibl,is) = p.cti;
            end
			p.thet(ibl,is) = thi;
			p.dstr(ibl,is) = dsi;
			p.uedg(ibl,is) = p.uei;
			p.mass(ibl,is) = dsi*p.uei;
			p.tau(ibl,is)  = 0.5*p.r2*p.u2*p.u2*p.cf2;
			p.dis(ibl,is)  =     p.r2*p.u2*p.u2*p.u2*p.di2*p.hs2*0.5;
			p.ctq(ibl,is)  = p.cq2;

			%//------ set "1" variables to "2" variables for next streamwise station
			p=blprv(p.xsi,p.ami,p.cti,thi,dsi,p.dswaki,p.uei,p);
			p=blkin(p);

			p=stepbl(p);

			%//------ turbulent intervals will follow transition interval or te
			if(p.tran || ibl==p.iblte(is)) 
				p.turb = true;
				
				%//------- save transition location
				p.tforce(is) = p.trforc;
				p.xssitr(is) = p.xt;
            end

			p.tran = false;
            break
                end
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
            end

			
% 			disp('     mrchdu: convergence failed a');
% 			if(m_bTrace)pXFile->WriteString(str);

			if (p.dmax<= 0.1) %goto stop109;	
                
                
                
                
                
                p=blprv(p.xsi,p.ami,p.cti,thi,dsi,p.dswaki,p.uei,p);
			p=blkin(p);
			
			%//------- check for transition and set appropriate flags and things
			if((p.simi==false) && (p.turb==false)) 
				p=trchek(p);
				p.ami = p.ampl2;
				if( p.tran) 
                    p.itran(is) = ibl;
                end
				if(p.tran==false)
                    p.itran(is) = ibl+2;
                end
            end
			
			%//------- set all other extrapolated values for current station
			if(ibl<p.itran(is))  
                p=blvar(1,p);
            end
			if(ibl>=p.itran(is)) 
                p=blvar(2,p);
            end
			if(p.wake) 
                p=blvar(3,p);
            end
			
			if(ibl<p.itran(is))
                p=blmid(1,p);
            end
			if(ibl>=p.itran(is))  
                p=blmid(2,p);
            end
			if(p.wake) 
                p=blmid(3,p);
            end
			
			%//------ pick up here after the newton iterations
% stop110:
			sens = p.sennew;

			%//------ store primary variables
			if(ibl<p.itran(is))
                p.ctau(ibl,is) = p.ami;
            else
                p.ctau(ibl,is) = p.cti;
            end
			p.thet(ibl,is) = thi;
			p.dstr(ibl,is) = dsi;
			p.uedg(ibl,is) = p.uei;
			p.mass(ibl,is) = dsi*p.uei;
			p.tau(ibl,is)  = 0.5*p.r2*p.u2*p.u2*p.cf2;
			p.dis(ibl,is)  =     p.r2*p.u2*p.u2*p.u2*p.di2*p.hs2*0.5;
			p.ctq(ibl,is)  = p.cq2;

			%//------ set "1" variables to "2" variables for next streamwise station
			p=blprv(p.xsi,p.ami,p.cti,thi,dsi,p.dswaki,p.uei,p);
			p=blkin(p);

			p=stepbl(p);

			%//------ turbulent intervals will follow transition interval or te
			if(p.tran || ibl==p.iblte(is)) 
				p.turb = true;
				
				%//------- save transition location
				p.tforce(is) = p.trforc;
				p.xssitr(is) = p.xt;
            end

			p.tran = false;
            continue
            end
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
				%//------ the current unconverged solution might still be reasonable...
				
			if(p.dmax > 0.1) 
				%//------- the current solution is garbage --> extrapolate values instead
				if(ibl>3) 
					if(ibl<=p.iblte(is)) 
						thi = p.thet(ibm,is) * sqrt(p.xssi(ibl,is)/p.xssi(ibm,is));
						dsi = p.dstr(ibm,is) * sqrt(p.xssi(ibl,is)/p.xssi(ibm,is));
						p.uei = p.uedg(ibm,is);
					
                    else
						if(ibl==p.iblte(is)+1) 
							p.cti = p.cte;
							thi = p.tte;
							dsi = p.dte;
							p.uei = p.uedg(ibm,is);
						
                        else
							thi = p.thet(ibm,is);
							ratlen = (p.xssi(ibl,is)-p.xssi(ibm,is)) / (10.0*p.dstr(ibm,is));
							dsi = (p.dstr(ibm,is) + thi*ratlen) / (1.0 + ratlen);
							p.uei = p.uedg(ibm,is);
                        end
                    end
					if(ibl==p.itran(is)) 
                        p.cti = 0.05;
                    end
					if(ibl>p.itran(is)) 
                        p.cti = p.ctau(ibm,is);
                    end
                end
				
            end

% stop109:
			p=blprv(p.xsi,p.ami,p.cti,thi,dsi,p.dswaki,p.uei,p);
			p=blkin(p);
			
			%//------- check for transition and set appropriate flags and things
			if((p.simi==false) && (p.turb==false)) 
				p=trchek(p);
				p.ami = p.ampl2;
				if( p.tran)
                    p.itran(is) = ibl;
                end
				if(p.tran==false) 
                    p.itran(is) = ibl+2;
                end
            end
			
			%//------- set all other extrapolated values for current station
			if(ibl<p.itran(is))  
                p=blvar(1,p);
            end
			if(ibl>=p.itran(is)) 
                p=blvar(2,p);
            end
			if(p.wake)  
                p=blvar(3,p);
            end
			if(ibl<p.itran(is)) 
                p=blmid(1,p);
            end
			if(ibl>=p.itran(is)) 
                p=blmid(2,p);
            end
			if(p.wake) 
                p=blmid(3,p);
            end
			%//------ pick up here after the newton iterations
% stop110:
			sens = p.sennew;

			%//------ store primary variables
			if(ibl<p.itran(is)) 
                p.ctau(ibl,is) = p.ami;
            else
                p.ctau(ibl,is) = p.cti;
            end
			p.thet(ibl,is) = thi;
			p.dstr(ibl,is) = dsi;
			p.uedg(ibl,is) = p.uei;
			p.mass(ibl,is) = dsi*p.uei;
			p.tau(ibl,is)  = 0.5*p.r2*p.u2*p.u2*p.cf2;
			p.dis(ibl,is)  =     p.r2*p.u2*p.u2*p.u2*p.di2*p.hs2*0.5;
			p.ctq(ibl,is)  = p.cq2;

			%//------ set "1" variables to "2" variables for next streamwise station
			p=blprv(p.xsi,p.ami,p.cti,thi,dsi,p.dswaki,p.uei,p);
			p=blkin(p);

			p=stepbl(p);

			%//------ turbulent intervals will follow transition interval or te
			if(p.tran || ibl==p.iblte(is)) 
				p.turb = true;
				
				%//------- save transition location
				p.tforce(is) = p.trforc;
				p.xssitr(is) = p.xt;
            end

			p.tran = false;
        end%//1000 continue
    end%// 2000 continue
% 	return true;
% }
end
function p=cdcalc(p)
% ////--------------------------------------
% //     cdcalc function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% {
% 	double p.thwake, p.urat, p.uewake, p.shwake, dx;
% 	int i,im;
	 sa = sin(p.alfa);
	 ca = cos(p.alfa);
	
	if(p.lvisc && p.lblini) 
		
% 		//---- set variables at the end of the wake
		p.thwake = p.thet(p.nbl(2),2);
		p.urat   = p.uedg(p.nbl(2),2)/p.qinf;
		p.uewake = p.uedg(p.nbl(2),2) * (1.0-p.tklam) / (1.0 - p.tklam*p.urat*p.urat);
		p.shwake = p.dstr(p.nbl(2),2)/p.thet(p.nbl(2),2);
		
% 		//---- extrapolate wake to downstream infinity using squire-young relation
% 		//      (reduces errors of the wake not being long enough)
		p.cd = 2.0*p.thwake * ((p.uewake/p.qinf)^(0.5*(5.0+p.shwake)));
	
    else
		p.cd = 0.0;
    end
	
% 	//--- calculate friction drag coefficient
	p.cdf = 0.0;
	for (is=1:2)
		for(ibl=3:p.iblte(is))
			i  = p.ipan(ibl,is);
			im = p.ipan(ibl-1,is);
			dx = (p.x(i) - p.x(im))*ca + (p.y(i) - p.y(im))*sa;
			p.cdf = p.cdf + 0.5*(p.tau(ibl,is)+p.tau(ibl-1,is))*dx * 2.0/p.qinf/p.qinf;
        end
    end
	
% 	return true;
end
function p=stmove(p)
% ////--------------------------------------
% //     stmove function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------


	%//--------------------------------------------------
	%//    moves stagnation point location to new panel.
	%//---------------------------------------------------
% 	int ibl, idif;
% 	double p.dudx;
	%//-- locate new stagnation point arc length sst from gam distribution
	p.istold = p.ist;
	p=stfind(p);
	
	if(p.istold==p.ist) 
		
		%//--- recalculate new arc length array
		p=xicalc(p);
	
	
    else
		
		%//       write(*,*) 'stmove: resetting stagnation point'
		
		%//--- set new bl position -> panel position  pointers
		p=iblpan(p);
		
		%//--- set new inviscid bl edge velocity uinv from qinv
		p=uicalc(p);
		
		%//--- recalculate new arc length array
		p=xicalc(p);
		
		%//--- set  bl position -> system line  pointers
		p=iblsys(p);
		
		if(p.ist>p.istold) 
			%//---- increase in number of points on top side (is=1)
			idif = p.ist-p.istold;
			
			
			p.itran(1) = p.itran(1) + idif;
			p.itran(2) = p.itran(2) - idif;
			
			%//---- move top side bl variables downstream
			for (ibl=p.nbl(1):-1:idif+2)
				
				p.ctau(ibl,1) = p.ctau(ibl-idif,1);
				p.thet(ibl,1) = p.thet(ibl-idif,1);
				p.dstr(ibl,1) = p.dstr(ibl-idif,1);
				p.uedg(ibl,1) = p.uedg(ibl-idif,1);
            end
			
			%//---- set bl variables between old and new stagnation point
			p.dudx = p.uedg(idif+2,1)/p.xssi(idif+2,1);
			for (ibl=idif+1:-1: 2)
			
				p.ctau(ibl,1) = p.ctau(idif+2,1);
				p.thet(ibl,1) = p.thet(idif+2,1);
				p.dstr(ibl,1) = p.dstr(idif+2,1);
				p.uedg(ibl,1) = p.dudx * p.xssi(ibl,1);
            end
			
			%//---- move bottom side bl variables upstream
			for (ibl=2: p.nbl(2))
				p.ctau(ibl,2) = p.ctau(ibl+idif,2);
				p.thet(ibl,2) = p.thet(ibl+idif,2);
				p.dstr(ibl,2) = p.dstr(ibl+idif,2);
				p.uedg(ibl,2) = p.uedg(ibl+idif,2);
            end
		
        else
			%//---- increase in number of points on bottom side (is=2)
			idif = p.istold-p.ist;
			
			p.itran(1) = p.itran(1) - idif;
			p.itran(2) = p.itran(2) + idif;
			
			%//---- move bottom side bl variables downstream
			for (ibl=p.nbl(2):-1:idif+2)
				
				p.ctau(ibl,2) = p.ctau(ibl-idif,2);
				p.thet(ibl,2) = p.thet(ibl-idif,2);
				p.dstr(ibl,2) = p.dstr(ibl-idif,2);
				p.uedg(ibl,2) = p.uedg(ibl-idif,2);
            end
			
			%//---- set bl variables between old and new stagnation point
			p.dudx = p.uedg(idif+2,2)/p.xssi(idif+2,2);
			for (ibl=idif+1:-1: 2)
				
				p.ctau(ibl,2) = p.ctau(idif+2,2);
				p.thet(ibl,2) = p.thet(idif+2,2);
				p.dstr(ibl,2) = p.dstr(idif+2,2);
				p.uedg(ibl,2) = p.dudx * p.xssi(ibl,2);
            end
			
			%//---- move top side bl variables upstream
			for(ibl=2: p.nbl(1))
				p.ctau(ibl,1) = p.ctau(ibl+idif,1);
				p.thet(ibl,1) = p.thet(ibl+idif,1);
				p.dstr(ibl,1) = p.dstr(ibl+idif,1);
				p.uedg(ibl,1) = p.uedg(ibl+idif,1);
            end
        end
		
    end
	
	%//-- set new p.mass array since ue has been tweaked
	for (is=1: 2)
		for(ibl=2:p.nbl(is)) 
			p.mass(ibl,is) = p.dstr(ibl,is)*p.uedg(ibl,is);
        end
    end
	
% 	return true; 
% }
end
function p=qvfue(p)
% ////--------------------------------------
% //     qvfue function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% 	//-------------------------------------------------------------
% 	//    sets panel viscous tangential velocity from viscous ue
% 	//-------------------------------------------------------------
% 	
	for ( is=1:2)
		for ( ibl=2:p.nbl(is))
			 i = p.ipan(ibl,is);
			p.qvis(i) = p.vti(ibl,is)*p.uedg(ibl,is);
        end
    end
	
% 	return true;
end
function p=gamqv(p)
% ////--------------------------------------
% //     gamqv function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------


	for ( i=1:p.n)
        p.gam(i)   = p.qvis(i);
        p.gam_a(i) = p.qinv_a(i);
    end
	
% 	return true;
% }
end
function p=update(p)
% ////--------------------------------------
% //     update function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

%//------------------------------------------------------------------
%//      adds on newton deltas to boundary layer variables.
%//      checks for excessive changes and underrelaxes if necessary.
%//      calculates max and rms changes.
%//      also calculates the change in the global variable "ac".
%//        if p.lalfa=true , "ac" is p.cl
%//        if p.lalfa=false, "ac" is alpha
%//------------------------------------------------------------------
		
% 	int i, ip, is, iv, iw, j, js, jv, ibl, jbl, kbl;
% 	double p.unew(IVX,3), p.u_ac(IVX,3);
% 	double p.qnew(IQX),p.q_ac(IQX);
% 	double p.dalmax, p.dalmin, p.dclmax, p.dclmin ,p.dx, p.dx_a, p.ag, p.ag_ms, p.ag_ac;
% 	double p.dac, p.dhi, p.dlo, p.dctau, p.dthet, p.dmass, p.duedg, p.ddstr;
% 	double p.dn1, p.dn2, p.dn3, p.dn4, rdn1,rdn2,rdn3,rdn4;
% 	double p.dswaki, p.hklim, p.msq, p.dsw;
% 	double p.dui, p.dui_ac, p.ue_m, p.uinv_ac,sa,ca, p.beta, p.beta_msq, p.bfac, p.bfac_msq;
% 	double p.clnew, p.cl_a, p.cl_ms, p.cl_ac, p.cginc;
% 	double p.cpg1,p.cpg1_ms, p.cpi_q, p.cpc_cpi, p.cpg1_ac,p.cpg2, p.cpg2_ms, p.cpg2_ac;
% 	CString p.vmxbl;
	
	%//---- max allowable alpha changes per iteration
	p.dalmax =  0.5*p.dtor;
	p.dalmin = -0.5*p.dtor;
%//---- max allowable p.cl change per iteration
	p.dclmax =  0.5;
	p.dclmin = -0.5;
	if(p.matyp~=1) 
        p.dclmin = max(-0.5, -0.9*p.cl) ;
    end
	p.hstinv = p.gamm1*(p.minf/p.qinf)*(p.minf/p.qinf) / (1.0 + 0.5*p.gamm1*p.minf*p.minf);

%//--- calculate new ue distribution assuming no under-relaxation
%//--- also set the sensitivity of ue wrt to alpha or re
	for (is=1: 2)
		for(ibl=2:p.nbl(is))
			i = p.ipan(ibl,is);
			p.dui    = 0.0;
			p.dui_ac = 0.0;
			for (js=1:2)
				for (jbl=2:p.nbl(js))
					j  = p.ipan(jbl,js);
					jv = p.isys(jbl,js);
					p.ue_m = -p.vti(ibl,is)*p.vti(jbl,js)*p.dij(i,j);
					p.dui    = p.dui    + p.ue_m*(p.mass(jbl,js)+p.vdel(3,1,jv));
					p.dui_ac = p.dui_ac + p.ue_m*(             -p.vdel(3,2,jv));
                end
            end
			
			%//------- p.uinv depends on "ac" only if "ac" is alpha
			if(p.lalfa) 
                p.uinv_ac = 0.0;
            else
                p.uinv_ac = p.uinv_a(ibl,is);
            end
			
			p.unew(ibl,is) = p.uinv(ibl,is) + p.dui;
			p.u_ac(ibl,is) = p.uinv_ac      + p.dui_ac;
			
        end
    end

	%//--- set new qtan from new ue with appropriate sign change

	for (is=1: 2)
		for(ibl=2:p.iblte(is))
			i = p.ipan(ibl,is);
			p.qnew(i) = p.vti(ibl,is)*p.unew(ibl,is);
			p.q_ac(i) = p.vti(ibl,is)*p.u_ac(ibl,is);
        end
    end
	
	%//--- calculate new p.cl from this new qtan
	sa = sin(p.alfa);
	ca = cos(p.alfa);

	
	p.beta = sqrt(1.0 - p.minf*p.minf);
	p.beta_msq = -0.5/p.beta;
	
	p.bfac     = 0.5*p.minf*p.minf / (1.0 + p.beta);
	p.bfac_msq = 0.5         / (1.0 + p.beta)...
			  - p.bfac        / (1.0 + p.beta) * p.beta_msq;
	
	p.clnew = 0.0;
	p.cl_a  = 0.0;
	p.cl_ms = 0.0;
	p.cl_ac = 0.0;
	
	i = 1;
	p.cginc = 1.0 - (p.qnew(i)/p.qinf)*(p.qnew(i)/p.qinf);
	p.cpg1  = p.cginc / (p.beta + p.bfac*p.cginc);
	p.cpg1_ms = -p.cpg1/(p.beta + p.bfac*p.cginc)*(p.beta_msq + p.bfac_msq*p.cginc);
	
	p.cpi_q = -2.0*p.qnew(i)/p.qinf/p.qinf;
	p.cpc_cpi = (1.0 - p.bfac*p.cpg1)/ (p.beta + p.bfac*p.cginc);
	p.cpg1_ac = p.cpc_cpi*p.cpi_q*p.q_ac(i);

	
	for (i=1:p.n)
		ip = i+1;
		if(i==p.n)
            ip = 1;
        end
		
		p.cginc = 1.0 - (p.qnew(ip)/p.qinf)*(p.qnew(ip)/p.qinf);
		p.cpg2  = p.cginc / (p.beta + p.bfac*p.cginc);
		p.cpg2_ms = -p.cpg2/(p.beta + p.bfac*p.cginc)*(p.beta_msq + p.bfac_msq*p.cginc);
		
		p.cpi_q = -2.0*p.qnew(ip)/p.qinf/p.qinf;
		p.cpc_cpi = (1.0 - p.bfac*p.cpg2)/ (p.beta + p.bfac*p.cginc);
		p.cpg2_ac = p.cpc_cpi*p.cpi_q*p.q_ac(ip);
		
		p.dx   =  (p.x(ip) - p.x(i))*ca + (p.y(ip) - p.y(i))*sa;
		p.dx_a = -(p.x(ip) - p.x(i))*sa + (p.y(ip) - p.y(i))*ca;
		
		p.ag    = 0.5*(p.cpg2    + p.cpg1   );
		p.ag_ms = 0.5*(p.cpg2_ms + p.cpg1_ms);
		p.ag_ac = 0.5*(p.cpg2_ac + p.cpg1_ac);
		
		p.clnew = p.clnew + p.dx  *p.ag;
		p.cl_a  = p.cl_a  + p.dx_a*p.ag;
		p.cl_ms = p.cl_ms + p.dx  *p.ag_ms;
		p.cl_ac = p.cl_ac + p.dx  *p.ag_ac;
		
		p.cpg1    = p.cpg2;
		p.cpg1_ms = p.cpg2_ms;
		p.cpg1_ac = p.cpg2_ac;
    end
	
	%//--- initialize under-relaxation factor
	rlx = 1.0;
	
	if(p.lalfa) 
		%//===== alpha is prescribed: ac is p.cl
		
		%//---- set change in re to account for p.cl changing, since re = re(p.cl)
		p.dac = (p.clnew - p.cl) / (1.0 - p.cl_ac - p.cl_ms*2.0*p.minf*p.minf_cl);
		
		%//---- set under-relaxation factor if re change is too large
		if(rlx*p.dac > p.dclmax) 
            rlx = p.dclmax/p.dac;
        end
		if(rlx*p.dac < p.dclmin) 
            rlx = p.dclmin/p.dac;
        end
	
    else
		%//===== p.cl is prescribed: ac is alpha
		
		%//---- set change in alpha to drive p.cl to prescribed value
		p.dac = (p.clnew - p.clspec) / (0.0 - p.cl_ac - p.cl_a);
		
		%//---- set under-relaxation factor if alpha change is too large
		if(rlx*p.dac > p.dalmax) 
            rlx = p.dalmax/p.dac;
        end
		if(rlx*p.dac < p.dalmin) 
            rlx = p.dalmin/p.dac;
        end
    end
	p.rmsbl = 0.0;
	p.rmxbl = 0.0;
	p.dhi = 1.5;
	p.dlo = -.5;
	%//--- calculate changes in bl variables and under-relaxation if needed

	for(is=1:2)
		for(ibl=2:p.nbl(is))
			iv = p.isys(ibl,is);
			%//------- set changes without underrelaxation
			p.dctau = p.vdel(1,1,iv) - p.dac*p.vdel(1,2,iv);
			p.dthet = p.vdel(2,1,iv) - p.dac*p.vdel(2,2,iv);
			p.dmass = p.vdel(3,1,iv) - p.dac*p.vdel(3,2,iv);
			p.duedg = p.unew(ibl,is) + p.dac*p.u_ac(ibl,is)  -  p.uedg(ibl,is);
			p.ddstr = (p.dmass - p.dstr(ibl,is)*p.duedg)/p.uedg(ibl,is);
			%//------- normalize changes
			if(ibl<p.itran(is)) 
                p.dn1 = p.dctau / 10.0;
            else
                p.dn1 = p.dctau / p.ctau(ibl,is);
            end
			p.dn2 = p.dthet / p.thet(ibl,is);
			p.dn3 = p.ddstr / p.dstr(ibl,is);
			p.dn4 = abs(p.duedg)/0.25;
			%//------- accumulate for rms change
			p.rmsbl = p.rmsbl + p.dn1*p.dn1 + p.dn2*p.dn2 + p.dn3*p.dn3 + p.dn4*p.dn4;
			%//------- see if p.ctau needs underrelaxation
			rdn1 = rlx*p.dn1;
			if(abs(p.dn1) > abs(p.rmxbl)) 
				p.rmxbl = p.dn1;
				if(ibl<p.itran(is)) 
                    p.vmxbl = 'n';
                end
				if(ibl>=p.itran(is)) 
                    p.vmxbl = 'c';
                end
				p.imxbl = ibl;
				p.ismxbl = is;
            end
			if(rdn1 > p.dhi)
                rlx = p.dhi/p.dn1;
            end
			if(rdn1 < p.dlo)
                rlx = p.dlo/p.dn1;
            end
			%//------- see if theta needs underrelaxation
			rdn2 = rlx*p.dn2;
			if(abs(p.dn2) > abs(p.rmxbl)) 
				p.rmxbl = p.dn2;
				p.vmxbl = 't';
				p.imxbl = ibl;
				p.ismxbl = is;
            end
			if(rdn2 > p.dhi)
                rlx = p.dhi/p.dn2;
            end
			if(rdn2 < p.dlo) 
                rlx = p.dlo/p.dn2;
            end
			%//------- see if dstar needs underrelaxation
			rdn3 = rlx*p.dn3;
			if(abs(p.dn3) > abs(p.rmxbl)) 
				p.rmxbl = p.dn3;
				p.vmxbl = 'd';
				p.imxbl = ibl;
				p.ismxbl = is;
            end
			if(rdn3 > p.dhi)
                rlx = p.dhi/p.dn3;
            end
			if(rdn3 < p.dlo) 
                rlx = p.dlo/p.dn3;
            end
			
			%//------- see if ue needs underrelaxation
			rdn4 = rlx*p.dn4;
			if(abs(p.dn4) > abs(p.rmxbl)) 
				p.rmxbl = p.duedg;
				p.vmxbl = 'u';
				p.imxbl = ibl;
				p.ismxbl = is;
            end
			if(rdn4 > p.dhi) 
                rlx = p.dhi/p.dn4;
            end
			if(rdn4 < p.dlo) 
                rlx = p.dlo/p.dn4;
            end
        end
    end

	
	%//--- set true rms change
	p.rmsbl = sqrt( p.rmsbl / (4.0.*( p.nbl(1)+p.nbl(2) )) );
	
	if(p.lalfa) 
		%//---- set underrelaxed change in reynolds number from change in lift
		p.cl = p.cl + rlx*p.dac;
	
    else
		%//---- set underrelaxed change in alpha
		p.alfa = p.alfa + rlx*p.dac;
		p.adeg = p.alfa/dtor;
    end
	
	%//--- update bl variables with underrelaxed changes
	for(is=1:2)
		for(ibl=2:p.nbl(is))
			iv = p.isys(ibl,is);
			
			p.dctau = p.vdel(1,1,iv) - p.dac*p.vdel(1,2,iv);
			p.dthet = p.vdel(2,1,iv) - p.dac*p.vdel(2,2,iv);
			p.dmass = p.vdel(3,1,iv) - p.dac*p.vdel(3,2,iv);
			p.duedg = p.unew(ibl,is) + p.dac*p.u_ac(ibl,is)  -  p.uedg(ibl,is);
			p.ddstr = (p.dmass - p.dstr(ibl,is)*p.duedg)/p.uedg(ibl,is);
			
			p.ctau(ibl,is) = p.ctau(ibl,is) + rlx*p.dctau;
			p.thet(ibl,is) = p.thet(ibl,is) + rlx*p.dthet;
			p.dstr(ibl,is) = p.dstr(ibl,is) + rlx*p.ddstr;
			p.uedg(ibl,is) = p.uedg(ibl,is) + rlx*p.duedg;
			
			if(ibl>p.iblte(is)) 
				iw = ibl - p.iblte(is);
				p.dswaki = p.wgap(iw);
			
            else
                p.dswaki = 0.0;
            end
			%//------- eliminate absurd transients
			if(ibl>=p.itran(is)) 
                p.ctau(ibl,is) = min(p.ctau(ibl,is), 0.25);
            end
			
			if(ibl<=p.iblte(is)) 
                p.hklim = 1.02;
            else
                p.hklim = 1.00005;
            end
			
			p.msq = p.uedg(ibl,is)*p.uedg(ibl,is)*p.hstinv...
				/ (p.gamm1*(1.0 - 0.5*p.uedg(ibl,is)*p.uedg(ibl,is)*p.hstinv));
			p.dsw = p.dstr(ibl,is) - p.dswaki;
			p.dsw=dslim(p.dsw,p.thet(ibl,is),p.msq,p.hklim);
			p.dstr(ibl,is) = p.dsw + p.dswaki;
			
			%//------- set new p.mass defect (nonlinear update)
			p.mass(ibl,is) = p.dstr(ibl,is) * p.uedg(ibl,is);
        end
    end

	
	%//--- equate upper wake arrays to lower wake arrays
	for(kbl=1:p.nbl(2)-p.iblte(2))
		p.ctau(p.iblte(1)+kbl,1) = p.ctau(p.iblte(2)+kbl,2);
		p.thet(p.iblte(1)+kbl,1) = p.thet(p.iblte(2)+kbl,2);
		p.dstr(p.iblte(1)+kbl,1) = p.dstr(p.iblte(2)+kbl,2);
		p.uedg(p.iblte(1)+kbl,1) = p.uedg(p.iblte(2)+kbl,2);
		p.tau(p.iblte(1)+kbl,1) =  p.tau(p.iblte(2)+kbl,2);
		p.dis(p.iblte(1)+kbl,1) =  p.dis(p.iblte(2)+kbl,2);
		p.ctq(p.iblte(1)+kbl,1) =  p.ctq(p.iblte(2)+kbl,2);
    end

	%//      equivalence (p.va(1,1,1), p.unew(1,1)) , (p.vb(1,1,1), p.qnew(1)  )
	%//      equivalence (p.va(1,1,IVX), p.u_ac(1,1)) , (p.vb(1,1,ivx), p.q_ac(1)  )
    IVX=p.n;
    
	for (kk = 1:p.n) 
		p.vb(kk,1,1)   = p.qnew(kk);
		p.vb(kk,1,end) = p.q_ac(kk);
    end

	for (is=1:2)
		for(ibl=2:p.nbl(is))
			p.va(ibl,is,1)   = p.unew(ibl,is);
			p.va(ibl,is,end) = p.u_ac(ibl,is);
        end
    end
% 	return true;
end
function p=blsolve(p)
% ////--------------------------------------
% //     blsolve function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------


	%//-----------------------------------------------------------------
	%//     custom solver for coupled viscous-inviscid newton system:
	%//
	%//       a  |  |  .  |  |  .  |    d       r       s
	%//       b  a  |  .  |  |  .  |    d       r       s
	%//       |  b  a  .  |  |  .  |    d       r       s
	%//       .  .  .  .  |  |  .  |    d   =   r - dre s
	%//       |  |  |  b  a  |  .  |    d       r       s
	%//       |  z  |  |  b  a  .  |    d       r       s
	%//       .  .  .  .  .  .  .  |    d       r       s
	%//       |  |  |  |  |  |  b  a    d       r       s
	%//
	%//      a, b, z  3x3  blocks containing linearized bl equation coefficients
	%//      |        3x1  vectors containing mass defect influence 
	%//                    coefficients on ue
	%//      d        3x1  unknown vectors (newton deltas for ctau, theta, m)
	%//      r        3x1  residual vectors
	%//      s        3x1  re influence vectors
	%//-----------------------------------------------------------------
	%//
% 	int iv, kv, ivp, k, l;

	ivte1 = p.isys(p.iblte(1),1);
	%//
	for (iv=1:p.nsys)
		%//
		ivp = iv + 1;
		%//
		%//====== invert va(iv) block
		%//
		%//------ normalize first row
		pivot = 1.0 / p.va(1,1,iv);
		p.va(1,2,iv) = p.va(1,2,iv) * pivot;
		for (l=iv:p.nsys) 
            p.vm(1,l,iv) = p.vm(1,l,iv)*pivot;
        end
		p.vdel(1,1,iv) = p.vdel(1,1,iv)*pivot;
		p.vdel(1,2,iv) = p.vdel(1,2,iv)*pivot;
		%//
		%//------ eliminate lower first column in va block
		for (k=2:3)
			vtmp = p.va(k,1,iv);
			p.va(k,2,iv) = p.va(k,2,iv) - vtmp*p.va(1,2,iv);
			for (l=iv:p.nsys) 
                p.vm(k,l,iv) = p.vm(k,l,iv) - vtmp*p.vm(1,l,iv);
            end
			p.vdel(k,1,iv) = p.vdel(k,1,iv) - vtmp*p.vdel(1,1,iv);
			p.vdel(k,2,iv) = p.vdel(k,2,iv) - vtmp*p.vdel(1,2,iv);
        end
		%//
		%//------ normalize second row
		pivot = 1.0 / p.va(2,2,iv);
		for (l=iv:p.nsys) 
            p.vm(2,l,iv) = p.vm(2,l,iv)*pivot;
        end
		p.vdel(2,1,iv) = p.vdel(2,1,iv)*pivot;
		p.vdel(2,2,iv) = p.vdel(2,2,iv)*pivot;
		%//
		%//------ eliminate lower second column in va block
		k = 3;
		vtmp = p.va(k,2,iv);
		for (l=iv:p.nsys) 
            p.vm(k,l,iv) = p.vm(k,l,iv) - vtmp*p.vm(2,l,iv);
        end
		p.vdel(k,1,iv) = p.vdel(k,1,iv) - vtmp*p.vdel(2,1,iv);
		p.vdel(k,2,iv) = p.vdel(k,2,iv) - vtmp*p.vdel(2,2,iv);
		
		%//------ normalize third row
		pivot = 1.0/p.vm(3,iv,iv);
		for (l=ivp:p.nsys) 
            p.vm(3,l,iv) = p.vm(3,l,iv)*pivot;
        end
		p.vdel(3,1,iv) = p.vdel(3,1,iv)*pivot;
		p.vdel(3,2,iv) = p.vdel(3,2,iv)*pivot;
		%//
		%//
		%//------ eliminate upper third column in va block
		 vtmp1 = p.vm(1,iv,iv);
		 vtmp2 = p.vm(2,iv,iv);
		for(l=ivp: p.nsys)
			p.vm(1,l,iv) = p.vm(1,l,iv) - vtmp1*p.vm(3,l,iv);
			p.vm(2,l,iv) = p.vm(2,l,iv) - vtmp2*p.vm(3,l,iv);
        end
		p.vdel(1,1,iv) = p.vdel(1,1,iv) - vtmp1*p.vdel(3,1,iv);
		p.vdel(2,1,iv) = p.vdel(2,1,iv) - vtmp2*p.vdel(3,1,iv);
		p.vdel(1,2,iv) = p.vdel(1,2,iv) - vtmp1*p.vdel(3,2,iv);
		p.vdel(2,2,iv) = p.vdel(2,2,iv) - vtmp2*p.vdel(3,2,iv);
		%//
		%//------ eliminate upper second column in va block
		vtmp = p.va(1,2,iv);
		for (l=ivp:p.nsys)
            p.vm(1,l,iv) = p.vm(1,l,iv) - vtmp*p.vm(2,l,iv);
        end
		p.vdel(1,1,iv) = p.vdel(1,1,iv) - vtmp*p.vdel(2,1,iv);
		p.vdel(1,2,iv) = p.vdel(1,2,iv) - vtmp*p.vdel(2,2,iv);
		%//
		%//
		if(iv~=p.nsys) 
			%//
			%//====== eliminate vb(iv+1) block, rows  1 -> 3
			for (k=1: 3)
				vtmp1 = p.vb(k, 1,ivp);
				vtmp2 = p.vb(k, 2,ivp);
				vtmp3 = p.vm(k,iv,ivp);
				for(l=ivp: p.nsys) 
                    p.vm(k,l,ivp) = p.vm(k,l,ivp)-(vtmp1*p.vm(1,l,iv)+ vtmp2*p.vm(2,l,iv)+vtmp3*p.vm(3,l,iv));
                end
				p.vdel(k,1,ivp) = p.vdel(k,1,ivp)-(vtmp1*p.vdel(1,1,iv)+vtmp2*p.vdel(2,1,iv)+ vtmp3*p.vdel(3,1,iv));
				p.vdel(k,2,ivp) = p.vdel(k,2,ivp)-(vtmp1*p.vdel(1,2,iv)+vtmp2*p.vdel(2,2,iv)+ vtmp3*p.vdel(3,2,iv));
            end
			%//
			if(iv==ivte1)
				%//------- eliminate vz block
				ivz = p.isys(p.iblte(2)+1,2);
				%//
				for(k=1:3)
					vtmp1 = p.vz(k,1);
					vtmp2 = p.vz(k,2);
					for (l=ivp: p.nsys)
						p.vm(k,l,ivz) = p.vm(k,l,ivz)-(vtmp1*p.vm(1,l,iv)+ vtmp2*p.vm(2,l,iv));
                    end
					p.vdel(k,1,ivz) = p.vdel(k,1,ivz)-(vtmp1*p.vdel(1,1,iv)+ vtmp2*p.vdel(2,1,iv));
					p.vdel(k,2,ivz) = p.vdel(k,2,ivz)-(vtmp1*p.vdel(1,2,iv)+ vtmp2*p.vdel(2,2,iv));
                end
            end
                
			%//
			if(ivp~=p.nsys) 
				%//
				%//====== eliminate lower vm column
				for(kv=iv+2:p.nsys)
					vtmp1 = p.vm(1,iv,kv);
					vtmp2 = p.vm(2,iv,kv);
					vtmp3 = p.vm(3,iv,kv);
					%//
					if(abs(vtmp1)>p.vaccel)
						for(l=ivp:p.nsys) 
                            p.vm(1,l,kv) = p.vm(1,l,kv) - vtmp1*p.vm(3,l,iv);
                        end
						p.vdel(1,1,kv) = p.vdel(1,1,kv) - vtmp1*p.vdel(3,1,iv);
						p.vdel(1,2,kv) = p.vdel(1,2,kv) - vtmp1*p.vdel(3,2,iv);
                    end
					%//
					if(abs(vtmp2)>p.vaccel) 
						for (l=ivp:p.nsys) 
                            p.vm(2,l,kv) = p.vm(2,l,kv) - vtmp2*p.vm(3,l,iv);
                        end
						p.vdel(2,1,kv) = p.vdel(2,1,kv) - vtmp2*p.vdel(3,1,iv);
						p.vdel(2,2,kv) = p.vdel(2,2,kv) - vtmp2*p.vdel(3,2,iv);
                    end
					%//
					if(abs(vtmp3)>p.vaccel) 
						for(l=ivp:p.nsys) 
                            p.vm(3,l,kv) = p.vm(3,l,kv) - vtmp3*p.vm(3,l,iv);
                        end
						p.vdel(3,1,kv) = p.vdel(3,1,kv) - vtmp3*p.vdel(3,1,iv);
						p.vdel(3,2,kv) = p.vdel(3,2,kv) - vtmp3*p.vdel(3,2,iv);
                    end
					%//
                end
            end
        end
    end%//1000

	%//
	for (iv=p.nsys:-1:2)
		%//------ eliminate upper vm columns
		vtmp = p.vdel(3,1,iv);
		for (kv=iv-1:-1:1)
			p.vdel(1,1,kv) = p.vdel(1,1,kv) - p.vm(1,iv,kv)*vtmp;
			p.vdel(2,1,kv) = p.vdel(2,1,kv) - p.vm(2,iv,kv)*vtmp;
			p.vdel(3,1,kv) = p.vdel(3,1,kv) - p.vm(3,iv,kv)*vtmp;
        end
		vtmp = p.vdel(3,2,iv);
		for (kv=iv-1:-1:1)
			p.vdel(1,2,kv) = p.vdel(1,2,kv) - p.vm(1,iv,kv)*vtmp;
			p.vdel(2,2,kv) = p.vdel(2,2,kv) - p.vm(2,iv,kv)*vtmp;
			p.vdel(3,2,kv) = p.vdel(3,2,kv) - p.vm(3,iv,kv)*vtmp;
        end
		%//
    end
% 	return true;
end
function [p, m_cls,r_cls]=mrcl(cls,p)
% ////--------------------------------------
% //     mrcl function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% %     //-------------------------------------------
%     //     sets actual mach, reynolds numbers
%     //     from unit-cl values and specified cls
%     //     depending on matyp,retyp flags.
%     //-------------------------------------------
% 	double rrat;
	cla = max(cls, 0.000001);
	if(p.retyp<1 || p.retyp>3) 
% //		AfxMessageBox("mrcl:  illegal Re(cls) dependence trigger, Setting fixed Re ", MB_ICONSTOP | MB_OK);
% 		CString str;
% 		str.Format("    mrcl:  illegal Re(cls) dependence trigger, Setting fixed Re ");
% 		if(m_bTrace)pXFile->WriteString(str);
		p.retyp = 1; 
    end
	if(p.matyp<1 || p.matyp>3) 
% //		AfxMessageBox("mrcl:  illegal Mach(cls) dependence trigger\n Setting fixed mach", MB_ICONSTOP | MB_OK);
% 		CString str;
% 		str.Format("    mrcl:  illegal Mach(cls) dependence trigger\r\n Setting fixed Mach");
% 		if(m_bTrace)pXFile->WriteString(str);
		p.matyp = 1;
    end
	
	switch(p.matyp)  
		case 1 
			p.minf  = p.minf1;
			m_cls = 0.0;
			
		   
		case 2      
			p.minf  =  p.minf1/sqrt(cla);
			m_cls = -0.5*p.minf/cla;
			
	   
		case 3
			p.minf  = p.minf1;
			m_cls = 0.0;
			
		
    end
	
	switch(p.retyp)
		case 1
			p.reinf = p.reinf1;
			r_cls = 0.0;
			
		
		case 2    
			p.reinf =  p.reinf1/sqrt(cla);
			r_cls = -0.5*p.reinf/cla;
			
		
		case 3
			p.reinf =  p.reinf1/cla;
			r_cls = -p.reinf /cla;
			
		
		
    end
	if(p.minf >= 0.99) 
% //		AfxMessageBox("mrcl: cl too low for chosen mach(cl) dependence", MB_ICONSTOP | MB_OK);
% 		//TRACE("      artificially limiting mach to  0.99\n");
% 		CString str;
% 		str.Format("mrcl: Cl too low for chosen Mach(Cl) dependence\r\n");
% 		if(m_bTrace)pXFile->WriteString(str);
% 		str.Format("      artificially limiting mach to  0.99");
% 		if(m_bTrace)pXFile->WriteString(str);
		p.minf = 0.99;
		m_cls = 0.0;
    end
	
	rrat = 1.0;
	if(p.reinf1 > 0.0) 
        rrat = p.reinf/p.reinf1;
    end
	if(rrat > 100.0)  
% //		AfxMessageBox("mrcl: cl too low for chosen Re(Cl) dependence ", MB_ICONSTOP | MB_OK);
% 		//TRACE("     artificially limiting re to %f\n",reinf1*100.0);
% 		CString str;
% 		str.Format("mrcl: cl too low for chosen Re(Cl) dependence\r\n");
% 		if(m_bTrace)pXFile->WriteString(str);
% 		str.Format("      artificially limiting Re to %.0f\r\n",reinf1*100.0);
% 		if(m_bTrace)pXFile->WriteString(str);
		p.reinf = p.reinf1*100.0;
		r_cls = 0.0;
    end
% 	return true;
end
function p=comset(p)
% ////--------------------------------------
% //     comset function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------


% 	//---- set karman-tsien parameter p.tklam
	 p.beta = sqrt(1.0 - p.minf*p.minf);
	 p.beta_msq = -0.5/p.beta;
	
	p.tklam   = p.minf*p.minf / (1.0 + p.beta)/ (1.0 + p.beta);
	p.tkl_msq = 	1.0 / (1.0 + p.beta)/ (1.0 + p.beta)...
		- 2.0*p.tklam/ (1.0 + p.beta) * p.beta_msq;
	
% 	//---- set sonic pressure coefficient and speed
	if(p.minf==0.0) 
		p.cpstar = -999.0;
		p.qstar = 999.0;
	
    else
		p.cpstar = 2.0 / (p.gamma*p.minf*p.minf) *...
				((((1.0 + 0.5*p.gamm1*p.minf*p.minf)/(1.0 + 0.5*p.gamm1))^(p.gamma/p.gamm1)) ...
					- 1.0 );
		p.qstar = p.qinf/p.minf * sqrt( (1.0 + 0.5*p.gamm1*p.minf*p.minf)...
			/(1.0 + 0.5*p.gamm1		) );
    end
	
% 	return true;
% }
end
function p=initialize()
% ////--------------------------------------
% //     initialize function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

%//---------------------------------------------------
%//     variable initialization/default routine.
%//---------------------------------------------------
        
%         p.pi = 4.0*atan(1.0);
        p.hopi = 0.50/pi;
        p.qopi = 0.25/pi;
        p.dtor = pi/180.0;

        p.n=0;%// arcds : so that current airfoil is not initialized

        %(aijpiv, 0, sizeof(aijpiv));
        %(apanel, 0, sizeof(apanel));
        %(blsav, 0, sizeof(blsav));
        %(aij, 0, sizeof(aij));
        %(bij, 0, sizeof(bij));
        %(cij, 0, sizeof(cij));
        %(cpi, 0, sizeof(cpi));
        %(cpv, 0, sizeof(cpv));
        %(ctau, 0, sizeof(ctau));
        %(ctq, 0, sizeof(ctq));
        %(dij, 0, sizeof(dij));
        %(dis, 0, sizeof(dis));
        %(dq, 0, sizeof(dq));
        %(dqdg, 0, sizeof(dqdg));
        %(dqdm, 0, sizeof(dqdm));	memset(delt, 0, sizeof(delt));
        %(dstr, 0, sizeof(dstr));
        %(dzdg, 0, sizeof(dzdg));
        %(dzdm, 0, sizeof(dzdm));
        %(dzdn, 0, sizeof(dzdn));
        %(guxd, 0, sizeof(guxd));
        %(guxq, 0, sizeof(guxq));
        %(iblte, 0, sizeof(iblte));
        %(ipan, 0, sizeof(ipan));
        %(isys, 0, sizeof(isys));
        %(itran, 0, sizeof(itran));
        %(mass, 0, sizeof(mass));
        %(nbl, 0, sizeof(nbl));
        %(nx, 0, sizeof(nx));
        %(ny, 0, sizeof(ny));
        %(gamu, 0, sizeof(gamu));
        %(gam, 0, sizeof(gam));
        %(gam_a, 0, sizeof(gam_a));
        %(q, 0, sizeof(q));
        %(qf0, 0, sizeof(qf0));
        %(qf1, 0, sizeof(qf1));
        %(qf2, 0, sizeof(qf2));
        %(qf3, 0, sizeof(qf3));
        %(qinv, 0, sizeof(qinv));
        %(qinvu, 0, sizeof(qinvu));
        %(qinv_a, 0, sizeof(qinv_a));
        %(qvis, 0, sizeof(qvis));
        %(s, 0, sizeof(x));
        %(sb, 0, sizeof(xb));
        %(sig, 0, sizeof(sig));
        %(snew, 0, sizeof(snew));
        %(sig, 0, sizeof(sig));
        %(tau, 0, sizeof(tau));
        %(thet, 0, sizeof(thet));
        %(uedg, 0, sizeof(uedg));
        %(uinv, 0, sizeof(uinv));
        %(uslp, 0, sizeof(uslp));
        %(vti, 0, sizeof(vti));
        %(x, 0, sizeof(x));
        %(xb, 0, sizeof(xb));
        %(xbp, 0, sizeof(xbp));
        %(xp, 0, sizeof(xp));
        %(xssi, 0, sizeof(xssi));
        %(y, 0, sizeof(y));
        %(yb, 0, sizeof(yb));
        %(ybp, 0, sizeof(ybp));
        %(yp, 0, sizeof(yp));
        %(wgap, 0, sizeof(wgap));
        %(va, 0, sizeof(va));
        %(vb, 0, sizeof(vb));
        %(vdel, 0, sizeof(vdel));
        %(vm, 0, sizeof(vm));
        %(vs1, 0, sizeof(vs1));
        %(vs2, 0, sizeof(vs2));
        %(vsrez, 0, sizeof(vsrez));
        %(vsr, 0, sizeof(vsr));
        %(vsm, 0, sizeof(vsm));
        %(vsx, 0, sizeof(vsx));
        %(vz, 0, sizeof(vz));
        %(w1, 0, sizeof(w1));
        %(w2, 0, sizeof(w2));
        %(w3, 0, sizeof(w3));
        %(w4, 0, sizeof(w4));
        %(w5, 0, sizeof(w5));
        %(w6, 0, sizeof(w6));
        %(w7, 0, sizeof(w7));
        %(w8, 0, sizeof(w8));

%         p.//mdes
        %(wc, 0, sizeof(wc));
        %(sc, 0, sizeof(sc));
        %(scold, 0, sizeof(scold));
        %(xcold, 0, sizeof(xcold));
        %(ycold, 0, sizeof(ycold));
        %(sspec, 0, sizeof(sspec));
        %(xspoc, 0, sizeof(xspoc));
        %(yspoc, 0, sizeof(yspoc));
        %(qgamm, 0, sizeof(qgamm));
        %(qspec, 0, sizeof(qspec));
        %(qspecp, 0, sizeof(qspecp));
        %(alqsp, 0, sizeof(alqsp));
        %(clqsp, 0, sizeof(clqsp));
        %(cmqsp, 0, sizeof(cmqsp));

        p.agte = 0.0;
        p.ag0 = 0.0;
        p.qim0 = 0.0;
        p.qimold = 0.0;
        p.ssple = 0.0;
        p.dwc = 0.0;
        p.algam = 0.0;
        p.clgam = 0.0;
        p.cmgam = 0.0;

        p.niterq = 6;

        %//---- default cp/cv (air)
        p.gamma = 1.4;
        p.gamm1 = p.gamma - 1.0;

        %//---- set unity freestream speed
        p.qinf = 1.0;

        p.psio = 0.0;

        p.cl = 0.0;
        p.cm = 0.0;
        p.cd = 0.0;

        p.sigte = 0.0;
        p.gamte = 0.0;
        %//	sigte_a = 0.0;
        %//	gamte_a = 0.0;

        p.nsp = 0;
        p.nqsp = 0;

        p.awake = 0.0;
        p.avisc = 0.0;

        %//	kimage = 1;
        p.yimage = -10.0;
        p.limage = false;

        p.liqset = false; %// ???
        p.lgamu  = false;
        p.lqinu  = false;%//???
        p.lvisc  = false;
        p.lwake  = false;
        %//	lpacc  = false;
        p.lblini = false;
        p.lipan  = false;
        p.lqaij  = false;
        p.ladij  = false;
        p.lwdij  = false;
        p.lcpxx  = false;
        %//	lqvdes = false;
        p.lqspec = false;
        %//	lqrefl = false;
        p.lvconv = false;
        %//	lcpref = false;
        %//	lforef = false;
        %//	lpfile = false;
        %//	lpfilx = false;
        %//	lppsho = false;
        p.leiw   = false;
        p.lscini = false;

        %//	lclip  = false;
        %//	lvlab  = true;
        %//	lcminp = false;
        %//	lhmomp = false;

        %//	lcurs  = true;
        %//	lland  = true;
        p.lgsame = false;

        %//	lgparm = true;
        %//	lplcam = false;

        p.sharp = false;
        p.lalfa = false;
        p.lbflap = false;
        p.lflap = false;
        p.trforc = false;
        p.simi = false;
        p.tran = false;
        p.turb = false;
        p.wake = false;
        p.trfree = false;
        p.tforce(0+1) =false;
        p.tforce(1+1) =false;
        p.tforce(2+1) =false;

        p.thickb = 0.0;
        p.cambrb = 0.0;

        %//---- input airfoil will not be normalized
        %//	lnorm = false;

        %//---- airfoil will not be forced symmetric
        p.lqsym = false;
        %//	lgsym = false;

        %//---- endpoint slopes will be matched
        p.lqslop = true;
        %//	lgslop = true;
        %//	lcslop = true;

        %//---- buffer and current airfoil flap hinge coordinates
        p.xbf = 0.0;
        p.ybf = 0.0;
        p.xof = 0.0;
        p.yof = 0.0;

        %//	ncpref = 0;
        %//        p.        p.        p.        p.       n
        %//---- circle plane array size (largest 2  + 1 that will fit array size)
        IQX=360;
        p.ann = log(((2*IQX)-1))/log(2.0);
        p.nn = ( p.ann + 0.00001 );
        p.tmp = 1;

        for (l=0:p.nn)
        p.tmp = 2*p.tmp;
        end
        p.nc1 = p.tmp + 1;
        %//	nc1 = (int)pow(2,nn) + 1;
        ICX=257;
        if(p.nc1 > ICX) 
            p.tmp = 1;
            for (l=0:p.nn-1)
                p.tmp = 2*p.tmp;
            end
            p.nc1 = p.tmp+1;
        %//		nc1 = pow(2,(nn-1)) + 1; //257 instead of ICX in original source code
        end
        p.IQX=360;
for I=1: p.IQX
p.gamu(I,1) = 0.;
p.gamu(I,2) = 0.;
p.gam(I) = 0.;
p.gam_a(I) = 0.;
end
p.IWX=p.IQX/8+2;
p.IZX=p.IQX+p.IWX;
for I=1: p.IZX
p.sig(I) = 0.;
end


        %//---- default cm reference location
        p.xcmref = 0.25;
        p.ycmref = 0.0;

        p.xoctr(1) = 1.0;
        p.xoctr(2) = 1.0;
        p.yoctr(1) = 0.0;
        p.yoctr(2) = 0.0;
        p.waklen = 1.0;

%         p.//added arcds : no wake yet
        p.nw = 0;

%         p.//added arcds : no flap yet
        p.hmom = 0.0;
        p.hfx  = 0.0;
        p.hfy  = 0.0;

%         p.//added arcds : fortran initializes to 0
        p.imxbl  = 0;
        p.ismxbl = 0;
        p.ist = 0;
        p.nb =0;


        p.dwte = 0.0;
        p.qinfbl = 0.0;
        p.tkbl = 0.0;
        p.tkbl_ms = 0.0;
        p.rstbl = 0.0;
        p.rstbl_ms = 0.0;
        p.hstinv = 0.0;
        p.hstinv_ms = 0.0;
        p.reybl = 0.0;
        p.reybl_ms = 0.0;
        p.reybl_re = 0.0;
        p.gambl = 0.0;
        p.gm1bl = 0.0;
        p.hvrat = 0.0;
        p.bule = 0.0;
        p.xiforc = 0.0;
        p.amcrit = 0.0;
        p.x2 = 0.0;
        p.u2 = 0.0;
        p.t2 = 0.0;
        p.d2 = 0.0;
        p.s2 = 0.0;
        p.ampl2 = 0.0;
        p.u2_uei = 0.0;
        p.u2_ms = 0.0;
        p.dw2 = 0.0;
        p.h2 = 0.0;
        p.h2_t2 = 0.0;
        p.h2_d2 = 0.0;
        p.m2 = 0.0;
        p.m2_u2 = 0.0;
        p.m2_ms = 0.0;
        p.r2 = 0.0; 
        p.r2_u2 = 0.0;
        p.r2_ms = 0.0;
        p.v2 = 0.0;
        p.v2_u2 = 0.0;
        p.v2_ms = 0.0;
        p.v2_re = 0.0;
        p.hk2 = 0.0;
        p.hk2_u2 = 0.0;
        p.hk2_t2 = 0.0;
        p.hk2_d2 = 0.0;
        p.hk2_ms = 0.0;
        p.hs2 = 0.0;
        p.hs2_u2 = 0.0;
        p.hs2_t2 = 0.0;
        p.hs2_d2 = 0.0;
        p.hs2_ms = 0.0;
        p.hs2_re = 0.0;
        p.hc2 = 0.0;
        p.hc2_u2 = 0.0;
        p.hc2_t2 = 0.0;
        p.hc2_d2 = 0.0;
        p.hc2_ms = 0.0;
        p.rt2 = 0.0;
        p.rt2_u2 = 0.0;
        p.rt2_t2 = 0.0;
        p.rt2_ms = 0.0;
        p.rt2_re = 0.0;
        p.cf2 = 0.0;
        p.cf2_u2 = 0.0;
        p.cf2_t2 = 0.0;
        p.cf2_d2 = 0.0;
        p.cf2_ms = 0.0;
        p.cf2_re = 0.0;
        p.di2 = 0.0;
        p.di2_u2 = 0.0;
        p.di2_t2 = 0.0;
        p.di2_d2 = 0.0;
        p.di2_s2 = 0.0;
        p.di2_ms = 0.0;
        p.di2_re = 0.0;
        p.us2 = 0.0;
        p.us2_u2 = 0.0;
        p.us2_t2 = 0.0;
        p.us2_d2 = 0.0;
        p.us2_ms = 0.0;
        p.us2_re = 0.0;
        p.cq2 = 0.0;
        p.cq2_u2 = 0.0;
        p.cq2_t2 = 0.0;
        p.cq2_d2 = 0.0;
        p.cq2_ms = 0.0;
        p.cq2_re = 0.0;
        p.de2 = 0.0;
        p.de2_u2 = 0.0;
        p.de2_t2 = 0.0;
        p.de2_d2 = 0.0;
        p.de2_ms = 0.0;
        p.x1 = 0.0;
        p.u1 = 0.0;
        p.t1 = 0.0;
        p.d1 = 0.0;
        p.s1 = 0.0;
        p.ampl1 = 0.0;
        p.u1_uei = 0.0;
        p.u1_ms = 0.0;
        p.dw1 = 0.0;
        p.h1 = 0.0;
        p.h1_t1 = 0.0;
        p.h1_d1 = 0.0;
        p.m1 = 0.0;
        p.m1_u1 = 0.0;
        p.m1_ms = 0.0;
        p.r1 = 0.0;
        p.r1_u1 = 0.0;
        p.r1_ms = 0.0;
        p.v1 = 0.0;
        p.v1_u1 = 0.0;
        p.v1_ms = 0.0;
        p.v1_re = 0.0;
        p.hk1 = 0.0;
        p.hk1_u1 = 0.0;
        p.hk1_t1 = 0.0;
        p.hk1_d1 = 0.0;
        p.hk1_ms = 0.0;
        p.hs1 = 0.0;
        p.hs1_u1 = 0.0;
        p.hs1_t1 = 0.0;
        p.hs1_d1 = 0.0;
        p.hs1_ms = 0.0;
        p.hs1_re = 0.0;
        p.hc1 = 0.0;
        p.hc1_u1 = 0.0;
        p.hc1_t1 = 0.0;
        p.hc1_d1 = 0.0;
        p.hc1_ms = 0.0;
        p.rt1 = 0.0;
        p.rt1_u1 = 0.0;
        p.rt1_t1 = 0.0;
        p.rt1_ms = 0.0;
        p.rt1_re = 0.0;
        p.cf1 = 0.0;
        p.cf1_u1 = 0.0;
        p.cf1_t1 = 0.0;
        p.cf1_d1 = 0.0;
        p.cf1_ms = 0.0;
        p.cf1_re = 0.0;
        p.di1 = 0.0;
        p.di1_u1 = 0.0;
        p.di1_t1 = 0.0;
        p.di1_d1 = 0.0;
        p.di1_s1 = 0.0;
        p.di1_ms = 0.0;
        p.di1_re = 0.0;
        p.us1 = 0.0;
        p.us1_u1 = 0.0;
        p.us1_t1 = 0.0;
        p.us1_d1 = 0.0;
        p.us1_ms = 0.0;
        p.us1_re = 0.0;
        p.cq1 = 0.0;
        p.cq1_u1 = 0.0;
        p.cq1_t1 = 0.0;
        p.cq1_d1 = 0.0;
        p.cq1_ms = 0.0;
        p.cq1_re = 0.0;
        p.de1 = 0.0;
        p.de1_u1 = 0.0;
        p.de1_t1 = 0.0;
        p.de1_d1 = 0.0;
        p.de1_ms = 0.0;
        p.xsf = 0.0;
        p.ysf = 0.0;
        p.cdp = 0.0;
        p.cdf = 0.0;
        p.alfa = 0.0;
        p.amax = 0.0;
        p.adeg = 0.0;
        p.rmxbl = 0.0;
        p.rmsbl = 0.0;
        p.rlx = 0.0;
        p.ante  = 0.0;
        p.ddef = 0.0;
        p.cpmn = 0.0;
        p.clspec = 0.0;
        p.minf = 0.0;
        p.reinf = 0.0;
        p.minf_cl = 0.0;
        p.reinf_cl = 0.0;

        p.sble = 0.0;
        p.chordb = 0.0;
        p.areab = 0.0;
        p.radble = 0.0;
        p.angbte = 0.0;
        p.ei11ba = 0.0;
        p.ei22ba = 0.0;
        p.apx1ba = 0.0;
        p.apx2ba = 0.0;
        p.ei11bt = 0.0;
        p.ei22bt = 0.0;
        p.apx1bt = 0.0;
        p.apx2bt = 0.0;
        p.sle = 0.0;
        p.xle = 0.0;
        p.yle = 0.0;
        p.xte = 0.0;
        p.yte = 0.0;
        p.chord = 0.0;
        p.ch = 0.0;
        p.cl_alf = 0.0;
        p.cl_msq = 0.0;
        p.cosa = 0.0;
        p.sina = 0.0;
        p.tklam = 0.0;
        p.tkl_msq = 0.0;
        p.cpstar = 0.0;
        p.qstar = 0.0;
        p.cpmni = 0.0;
        p.cpmnv = 0.0;
        p.xcpmni = 0.0;
        p.xcpmnv = 0.0;
        p.arad = 0.0;
        p.sst = 0.0;
        p.sst_go = 0.0;
        p.sst_gp = 0.0;
        p.dste = 0.0;
        p.aste = 0.0;

        p.qtan1 = 0.0;
        p.qtan2 = 0.0;
        p.z_qinf = 0.0;
        p.z_alfa = 0.0;
        p.z_qdof0 = 0.0;
        p.z_qdof1 = 0.0;
        p.z_qdof2 = 0.0;
        p.z_qdof3 = 0.0;
        p.cfm = 0.0;
        p.cfm_ms = 0.0;
        p.cfm_re = 0.0;
        p.cfm_u1 = 0.0;
        p.cfm_t1 = 0.0;
        p.cfm_d1 = 0.0;
        p.cfm_u2 = 0.0;
        p.cfm_t2 = 0.0;
        p.cfm_d2 = 0.0;
        p.xt = 0.0;
        p.xt_a1 = 0.0;
        p.xt_ms = 0.0;
        p.xt_re = 0.0;
        p.xt_xf = 0.0;
        p.xt_x1 = 0.0;
        p.xt_t1 = 0.0;
        p.xt_d1 = 0.0;
        p.xt_u1 = 0.0;
        p.xt_x2 = 0.0;
        p.xt_t2 = 0.0;
        p.xt_d2 = 0.0;
        p.xt_u2 = 0.0;

        p.npan = 140;
        p.cvpar = 1.0;
        p.cterat = 0.15;
        p.ctrrat = 0.2;

        %//---- default paneling refinement zone x/c endpoints
        p.xsref1 = 1.0;
        p.xsref2 = 1.0;
        p.xpref1 = 1.0;
        p.xpref2 = 1.0;

        %//---- drop tolerance for bl system solver
        p.vaccel = 0.01;
            p.m_bTrace = false;

    p.sccon = 5.6  ;
    p.gacon = 6.70;
    p.gbcon = 0.75;
    p.gbc0  = 0.60;
    p.gbc1  = 0.40;
    p.gccon = 18.0;
    p.dlcon =  0.9;
    p.ctcon = 0.01485111754659538130244;% //(ctcon = 0.5/(gacon**2 * gbcon))
    p.angtol = 40.0;

    %// fortran seems to initializes variables to 0
    p.mvisc = 0.0;

    %//initialize transition parameters until user changes them
    p.acrit     = 9.0;
    p.xstrip(1) = 1.0;
    p.xstrip(2) = 1.0;

    %//intialize analysis parameter  until user changes them
    %//---- default paneling parameters


    %//---- initialize freestream mach number to zero
    p.matyp = 1;
    p.minf1 = 0.0;

    %//---- drop tolerance for bl system solver

    %//---- default viscous parameters
    p.retyp = 1;
    p.reinf1 = 0.0;
end
function s=scalc(x,  y,    n)
% ////--------------------------------------
% //     scalc function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% //----------------------------------------
% //     calculates the arc length array s  |
% //     for a 2-d array of points (x,y).   |
% //----------------------------------------
	
	s(1) = 0.0;
	for (i=2:n)
		s(i) = s(i-1) + sqrt((x(i)-x(i-1))*(x(i)-x(i-1)) ...
								   +(y(i)-y(i-1))*(y(i)-y(i-1)));
    end
% 	return true;
% }
end
function xs=segspl( x,   s,  n)
% ////--------------------------------------
% //     segspl function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% //-----------------------------------------------
% //     splines x(s) array just like spline,      |
% //     but allows derivative discontinuities     |
% //     at segment joints.  segment joints are    |
% //     defined by identical successive s values. |
% //-----------------------------------------------

% 	 nseg, iseg0;
	
	if(s(1)==s(2)  ) 
        return ; %//stop 'segspl:  first input point duplicated'
    end
	if(s(n)==s(n-1)) 
        return ;% //stop 'segspl:  last  input point duplicated'
    end
	
	iseg0 = 1;
	for ( iseg=2:n-2)
		if(s(iseg)==s(iseg+1)) 
			nseg = iseg - iseg0 + 1;
% //			splind(x(iseg0),xs(iseg0),s(iseg0),nseg,-999.0,-999.0);
            x=[ones(iseg0-1,1) x];
            s=[ones(iseg0-1,1) s];
			[xs]=splind(x,s,nseg,-999.0,-999.0);
            xs=[ones(iseg0-1,1) xs];

			iseg0 = iseg+1;
        end
    end
	nseg = n - iseg0 + 1;
	
% //	splind(x(iseg0),xs(iseg0),s(iseg0),nseg,-999.0,-999.0);
            x=[ones(iseg0-1,1) x];
            s=[ones(iseg0-1,1) s];
	[xs]=splind(x,s,nseg,-999.0,-999.0);
                xs=[ones(iseg0-1,1) xs];

	
% 	return true;
end
function xs= splind( x,   s,  n,  xs1,  xs2)
% ////--------------------------------------
% //     splinf function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% 	 nmax=600;
% 	 a(601),b(601),c(601);
% 	//-------------------------------------------------------
% 	//     calculates spline coefficients for x(s).          |
% 	//     specified 1st derivative and/or usual zero 2nd    |
% 	//     derivative end conditions are used.               |
% 	
% 	//     to evaluate the spline at some value of s,        |
% 	//     use seval and/or deval.                           |
% 	//                                                       |
% 	//     s        independent variable array (input)       |
% 	//     x        dependent variable array   (input)       |
% 	//     xs       dx/ds array                (calculated)  |
% 	//     n        number of points           (input)       |
% 	//     xs1,xs2  endpoint derivatives       (input)       |
% 	//              if = 999.0, then usual zero second       |
% 	//              derivative end condition(s) are used     |
% 	//              if = -999.0, then zero third             |
% 	//              derivative end condition(s) are used     |
% 	//                                                       |
% 	//-------------------------------------------------------
% 
% 	 dsm, dsp;

% 	if(n>nmax) 
% //		AfxMessageBox("splind: array overflow, increase nmax", MB_ICONSTOP | MB_OK);
% 		CString str;
% 		str.Format("splind: array overflow, increase nmax");
% 		if(m_bTrace)pXFile->WriteString(str);
% 		return false;
% 	}
	for( i=2: n-1)
		dsm = s(i) - s(i-1);
		dsp = s(i+1) - s(i);
		b(i) = dsp;
		a(i) = 2.0*(dsm+dsp);
		c(i) = dsm;
		xs(i) = 3.0*((x(i+1)-x(i))*dsm/dsp + (x(i)-x(i-1))*dsp/dsm);
    end
	
	if(xs1>=998.0) 
% 		//----- set zero second derivative end condition
		a(1) = 2.0;
		c(1) = 1.0;
		xs(1) = 3.0*(x(2)-x(1)) / (s(2)-s(1));
    
	else 
		if(xs1<=-998.0) 
% 			//----- set zero third derivative end condition
			a(1) = 1.0;
			c(1) = 1.0;
			xs(1) = 2.0*(x(2)-x(1)) / (s(2)-s(1));
		
        else
% 			//----- set specified first derivative end condition
			a(1) = 1.0;
			c(1) = 0.0;
			xs(1) = xs1;
        end
    end
	
	
	if(xs2>=998.0) 
		b(n) = 1.0;
		a(n) = 2.0;
		xs(n) = 3.0*(x(n)-x(n-1)) / (s(n)-s(n-1));
	
    else
		if(xs2<=-998.0) 
			b(n) = 1.0;
			a(n) = 1.0;
			xs(n) = 2.0*(x(n)-x(n-1)) / (s(n)-s(n-1));
		
        else
			a(n) = 1.0;
			b(n) = 0.0;
			xs(n) = xs2;
        end
    end
	
	if(n==2 && xs1<=-998.0 && xs2<=-998.0) 
		b(n) = 1.0;
		a(n) = 2.0;
		xs(n) = 3.0*(x(n)-x(n-1)) / (s(n)-s(n-1));
    end
	
% 	//---- solve for derivative array xs
	xs=trisol(a,b,c,xs,n);
% 	return true;      
end
function d=trisol( a,  b,  c,  d,  kk)
% ////--------------------------------------
% //     trisol function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% 	//-----------------------------------------
% 	//     solves kk long, tri-diagonal system |
% 	//                                         |
% 	//             a c          d              |
% 	//             b a c        d              |
% 	//               b a .      .              |
% 	//                 . . c    .              |
% 	//                   b a    d              |
% 	//                                         |
% 	//     the righthand side d is replaced by |
% 	//     the solution.  a, c are destroyed.  |
% 	//-----------------------------------------
% 	 k;
	for (k=2:kk)
		 km = k-1;
		c(km) = c(km) / a(km);
		d(km) = d(km) / a(km);
		a(k) = a(k) - b(k)*c(km);
		d(k) = d(k) - b(k)*d(km);
    end
	
	d(kk) = d(kk)/a(kk);
	
	for(k=kk-1:-1: 1)
		d(k) = d(k) - c(k)*d(k+1);
    end
% 	return true;
% }
end
function [ t,  sle,   chord,area,  radle,  angte,ei11a,  ei22a,  apx1a,  apx2a,ei11t,  ei22t,  apx1t,  apx2t]=geopar( x,  xp,  y,  yp,  s,n)
% ////--------------------------------------
% //     geopar function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% 	 i;
% 	 chsq, curvle, ang1, ang2, xcena, ycena, slen, xcent, ycent;
% 	//------------------------------------------------------
% 	//     sets geometric parameters for airfoil shape
% 	//------------------------------------------------------
	sle=lefind(x,xp,y,yp,s,n);
	
	xle = seval(sle,x,xp,s,n);
	yle = seval(sle,y,yp,s,n);
	xte = 0.5*(x(1)+x(n));
	yte = 0.5*(y(1)+y(n));
	
	chsq = (xte-xle)*(xte-xle) + (yte-yle)*(yte-yle);
	chord = sqrt(chsq);
	
	curvle = curv(sle,x,xp,y,yp,s,n);
	
	radle = 0.0;
	if(abs(curvle) > 0.001*(s(n)-s(1))) 
        radle = 1.0 / curvle;
    end
	
	ang1 = atan2( -yp(1) , -xp(1) );
	ang2 = atanc(  yp(n) ,  xp(n) , ang1 );
	angte = ang2 - ang1;
	
	for (i=1:n)
        t(i) = 1.0;
    end
	
	[area,xcena,ycena,ei11a,ei22a,apx1a,apx2a]=aecalc(n,x,y,t, 1);
	[slen,xcent,ycent,ei11t,ei22t,apx1t,apx2t]=aecalc(n,x,y,t, 2);
	
% 	  thick, xthick, cambr, xcambr;
% 	 xcam(IQX), ycam(IQX), xthk(IQX), ythk(IQX), ycamp(IQX), ythkp(IQX);
% 	 ncam, nthk;
% 	//--- old, approximate thickness,camber routine (on discrete points only)
	[thick, xthick, cambr, xcambr]=tccalc(x,xp,y,yp,s,n );
% //--- more accurate thickness and camber estimates

% 	[xcam,ycam,ncam,xthk,ythk,nthk]=getcam(x,xp,y,yp,s,n );
% 	getmax(xcam,ycam,ycamp,ncam,xcambr,cambr);
% 	getmax(xthk,ythk,ythkp,nthk,xthick,thick);
% 	thick = 2.0*thick;

	thickb = thick;
	cambrb = cambr;

    
% 	//      write(*,1000) thick,xthick,cambr,xcambr
% 	// 1000 format( ' max thickness = ',f12.6,'  at x = ',f7.3,/' max camber    = ',f12.6,'  at x = ',f7.3)
	
% 	return true;
end	
function [area,xcen,  ycen,  ei11,  ei22,apx1,  apx2]=aecalc( n,  x,  y,  t,  itype)
% ////--------------------------------------
% //     aecalc function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------
  
% %//---------------------------------------------------------------
% %//     calculates geometric properties of shape x,y
% %//
% %//     input:
% %//       n      number of points
% %//       x(.)   shape coordinate point arrays
% %//       y(.)
% %//       t(.)   skin-thickness array, used only if itype = 2
% %//       itype  = 1 ...   integration is over whole area  dx dy
% %//              = 2 ...   integration is over skin  area   t ds
% %//
% %//     output:
% %//       xcen,ycen  centroid location
% %//       ei11,ei22  principal moments of inertia
% %//       apx1,apx2  principal-axis angles
% %//---------------------------------------------------------------
	
% 	 sint, aint, xint, yint, xxint, yyint, xyint;
% 	 eixx, eiyy, eixy, eisq;
% 	 dx, dy, xa, ya, ta, ds, da, c1, c2, sgn;
% 	 ip;
	sint  = 0.0;
	aint  = 0.0;
	xint  = 0.0;
	yint  = 0.0;
	xxint = 0.0;
	xyint = 0.0;
	yyint = 0.0;
	
	for ( io = 1:n)
		if(io==n) 
            ip = 1;
        
        else
            ip = io + 1;
        end
		
		
		dx =  x(io) - x(ip);
		dy =  y(io) - y(ip);
		xa = (x(io) + x(ip))*0.50;
		ya = (y(io) + y(ip))*0.50;
		ta = (t(io) + t(ip))*0.50;
		
		ds = sqrt(dx*dx + dy*dy);
		sint = sint + ds;
		
		if(itype==1) 
			%//-------- integrate over airfoil cross-section
			da = ya*dx;
			aint  = aint  +       da;
			xint  = xint  + xa   *da;
			yint  = yint  + ya   *da/2.0;
			xxint = xxint + xa*xa*da;
			xyint = xyint + xa*ya*da/2.0;
			yyint = yyint + ya*ya*da/3.0;
		
        else
			%//-------- integrate over skin thickness
			da = ta*ds;
			aint  = aint  +       da;
			xint  = xint  + xa   *da;
			yint  = yint  + ya   *da;
			xxint = xxint + xa*xa*da;
			xyint = xyint + xa*ya*da;
			yyint = yyint + ya*ya*da;
        end
    end
	
	area = aint;
	
	if(aint == 0.0) 
		xcen  = 0.0;
		ycen  = 0.0;
		ei11  = 0.0;
		ei22  = 0.0;
		apx1 = 0.0;
		apx2 = atan2(1.0,0.0);
% 		return false;
    end
	
	%//---- calculate centroid location
	xcen = xint/aint;
	ycen = yint/aint;
	
	%//---- calculate inertias
	eixx = yyint - (ycen)*(ycen)*aint;
	eixy = xyint - (xcen)*(ycen)*aint;
	eiyy = xxint - (xcen)*(xcen)*aint;
	
	%//---- set principal-axis inertias, ei11 is closest to "up-down" bending inertia
	eisq  = 0.25*(eixx - eiyy)*(eixx - eiyy)  + eixy*eixy;
	sgn = 1*sign( eiyy-eixx );
	ei11 = 0.5*(eixx + eiyy) - sgn*sqrt(eisq);
	ei22 = 0.5*(eixx + eiyy) + sgn*sqrt(eisq);
	
	if(ei11==0.0 || ei22==0.0) 
		%//----- vanishing section stiffness
		apx1 = 0.0;
		apx2 = atan2(1.0,0.0);
	
    else
		if(eisq/((ei11)*(ei22)) < ((0.001*sint)^4.0)) 
			%//----- rotationally-invariant section (circle, square, etc.)
			apx1 = 0.0;
			apx2 = atan2(1.0,0.0);
		
        else
			%//----- normal airfoil section
			c1 = eixy;
			s1 = eixx-ei11;
			
			c2 = eixy;
			s2 = eixx-ei22;
			
			if(abs(s1)>abs(s2)) 
				apx1 = atan2(s1,c1);
				apx2 = apx1 + 0.5*pi;
			
            else
				apx2 = atan2(s2,c2);
				apx1 = apx2 - 0.5*pi;
            end
			
			if(apx1<-0.5*pi) 
                apx1 = apx1 + pi;
            end
			if(apx1>+0.5*pi)
                apx1 = apx1 - pi;
            end
			if(apx2<-0.5*pi)
                apx2 = apx2 + pi;
            end
			if(apx2>+0.5*pi)
                apx2 = apx2 - pi;
            end
			
        end
    end
	
% 	return true;	
end
function sle=lefind( x,  xp,y,  yp,  s,  n)
% ////--------------------------------------
% //     lefind function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

%//------------------------------------------------------
%//     locates leading edge spline-parameter value sle
%//
%//     the defining condition is
%//         
%//      (x-xte,y-yte) . (x',y') = 0     at  s = sle
%//
%//     i.e. the surface tangent is normal to the chord
%//     line connecting x(sle),y(sle) and the te point.
%//------------------------------------------------------
% 	 i, iter;
% 	 dseps, dxte, dyte, dx, dy, dotp, dxds, dyds, dxdd, dydd;
% 	 res, ress, dsle;
% 	 xchord, ychord;
%//---- convergence tolerance
	dseps = (s(n)-s(1)) * 0.00001;
	
	%//---- set trailing edge point coordinates
	xte = 0.5*(x(1) + x(n));
	yte = 0.5*(y(1) + y(n));
	
	%//---- get first guess for sle
	for (i=3:n-2)
		dxte = x(i) - xte;
		dyte = y(i) - yte;
		dx = x(i+1) - x(i);
		dy = y(i+1) - y(i);
		dotp = dxte*dx + dyte*dy;
		if(dotp < 0.0) 
            break
        end
    end
	
	sle = s(i);
	
	%//---- check for sharp le case
	if(s(i) == s(i-1))
        return 
    end
	
	%//---- newton iteration to get exact sle value
	for (iter=1:50)
		xle  = seval(sle,x,xp,s,n);
		yle  = seval(sle,y,yp,s,n);
		dxds = deval(sle,x,xp,s,n);
		dyds = deval(sle,y,yp,s,n);
		dxdd = d2val(sle,x,xp,s,n);
		dydd = d2val(sle,y,yp,s,n);
		
		xchord = xle - xte;
		ychord = yle - yte;
		
		%//------ drive dot product between chord line and le tangent to zero
		res  = xchord*dxds + ychord*dyds;
		ress = dxds  *dxds + dyds  *dyds + xchord*dxdd + ychord*dydd;
		
		%//------ newton delta for sle 
		dsle = -res/ress;
		
		dsle = max( dsle , -0.02*abs(xchord+ychord) );
		dsle = min( dsle ,  0.02*abs(xchord+ychord) );
		sle = sle + dsle;
		if(abs(dsle) < dseps) 
            return 
        end
    end
	sle = s(i);
% 	return true;
end
function svl=seval( ss,  x,  xs,  s,  n)
% ////--------------------------------------
% //     seval function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% {
% 	//--------------------------------------------------
% 	//	  calculates x(ss)							   |
% 	//	   xs array must have been calculated by spline |
% 	//--------------------------------------------------
	 ilow = 1;
	 i = n;	
	
	while(i-ilow>1)
		 imid = round((i+ilow)/2);
		if(ss < s(imid)) 
            i = imid;
        else
            ilow = imid;
        end
    end
	
	 ds = s(i) - s(i-1);
	 t = (ss - s(i-1)) / ds;
	 cx1 = ds*xs(i-1) - x(i) + x(i-1);
	 cx2 = ds*xs(i)   - x(i) + x(i-1);
	 svl = t*x(i) + (1.0-t)*x(i-1) ...
					   + (t-t*t)*((1.0-t)*cx1 - t*cx2);
		
end
function deval=deval( ss,  x,  xs,  s,  n)
% ////--------------------------------------
% //     deval function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% 	//--------------------------------------------------
% 	//	   calculates dx/ds(ss) 						|
% 	//	   xs array must have been calculated by spline |
% 	//--------------------------------------------------
	 ilow = 1;
% //	 i = nc;
	 i = n; %///arcds modified
% 	 imid;
	
	
	while(i-ilow>1)
		imid = round((i+ilow)/2);
		if(ss < s(imid)) 
            i = imid;
        else
            ilow = imid;
        end
    end
	
	
	 ds = s(i) - s(i-1);
	 t = (ss - s(i-1)) / ds;
	 cx1 = ds*xs(i-1) - x(i) + x(i-1);
	 cx2 = ds*xs(i)	 - x(i) + x(i-1);
	 deval = x(i) - x(i-1) + (1.0-4.0*t+3.0*t*t)*cx1... 
								+ t*(3.0*t-2.0)*cx2;
	deval = deval/ds;
% 	return deval;
% }
end
function dtwoval=d2val( ss,  x,  xs,  s, n)
% ////--------------------------------------
% //     dtwoval function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------


% //--------------------------------------------------
% //     calculates d2x/ds2(ss)                       |
% //     xs array must have been calculated by spline |
% //--------------------------------------------------
% 	int imid;
	 ilow = 1;
	 i = n;
% stop10:
	while (i-ilow>1) 
        imid = round((i+ilow)/2);
	if(ss < s(imid))
        i = imid;
    else
        ilow = imid;
    end
    end
	
	

	 ds = s(i) - s(i-1);
	 t = (ss - s(i-1)) / ds;
	 cx1 = ds*xs(i-1) - x(i) + x(i-1);
	 cx2 = ds*xs(i)   - x(i) + x(i-1);
	 dtwoval = (6.0*t-4.0)*cx1 + (6.0*t-2.0)*cx2;
	dtwoval = dtwoval/ds/ds;
% 	return dtwoval;
end
function crv=curv( ss,  x,  xs,  y,  ys,  s,  n)
% ////--------------------------------------
% //     curv function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% //-----------------------------------------------
% //     calculates curvature of splined 2-d curve |
% //     at s = ss                                 |
% //                                               |
% //     s        arc length array of curve        |
% //     x, y     coordinate arrays of curve       |
% //     xs,ys    derivative arrays                |
% //              (calculated earlier by spline)   |
% //-----------------------------------------------
% 	 ilow, i, imid;
% 	 crv,ds, t, cx1, cx2, xd, xdd, cy1, cy2, yd, ydd, sd;
	
	
	ilow = 1;
	i = n;
	

	while(i-ilow>1)% goto stop11;
	imid = round((i+ilow)/2);
	if(ss < s(imid)) 
        i = imid;
    else
        ilow = imid;
    end
    end
	
% stop11:
	
	ds = s(i) - s(i-1);
	t = (ss - s(i-1)) / ds;
	
	cx1 = ds*xs(i-1) - x(i) + x(i-1);
	cx2 = ds*xs(i)   - x(i) + x(i-1);
	xd = x(i) - x(i-1) + (1.0-4.0*t+3.0*t*t)*cx1 + t*(3.0*t-2.0)*cx2;
	xdd = (6.0*t-4.0)*cx1 + (6.0*t-2.0)*cx2;
	
	cy1 = ds*ys(i-1) - y(i) + y(i-1);
	cy2 = ds*ys(i)   - y(i) + y(i-1);
	yd = y(i) - y(i-1) + (1.0-4.0*t+3.0*t*t)*cy1 + t*(3.0*t-2.0)*cy2;
	ydd = (6.0*t-4.0)*cy1 + (6.0*t-2.0)*cy2;
	
	sd = sqrt(xd*xd + yd*yd);
	sd = max(sd,0.001*ds);
	
	crv = (xd*ydd - yd*xdd) / sd/ sd/ sd;
	
% 	return crv;
end
function aa=atanc( y,  x,  thold)
% ////--------------------------------------
% //     atanc function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% //-------------------------------------------------------------
% //    atan2 function with branch cut checking.
% //
% //    increments position angle of point x,y from some previous
% //    value thold due to a change in position, ensuring that the
% //    position change does not cross the atan2 branch cut
% //    (which is in the -x direction).  for example:
% //
% //      atanc( -1.0 , -1.0 , 0.75*pi )  returns  1.25*pi , whereas
% //      atan2( -1.0 , -1.0 )            returns  -.75*pi .
% //
% //    typically, atanc is used to fill an array of angles:
% //
% //       theta(1) = atan2( y(1) , x(1) )
% //       do i=2, n
% //         theta[i] = atanc( y[i] , x[i] , theta(i-1) )
% //       end do
% //
% //    this will prevent the angle array theta(i) from jumping by 
% //    +/- 2 pi when the path x(i),y(i) crosses the negative x axis.
% //
% //    input:
% //      x,y     point position coordinates
% //      thold   position angle of nearby point
% //
% //    output:
% //      atanc   position angle of x,y
% //--------------------------------------------------------------
	
	 tpi = 2*pi;

% //---- set new position angle, ignoring branch cut in atan2 function for now
	
	 thnew = atan2( y , x );
	 dthet = thnew - thold;

% //---- angle change cannot exceed +/- pi, so get rid of any multiples of 2 pi 
	 dtcorr = dthet - tpi*round( (dthet + pi*sign(dthet))/tpi );
% 
% //---- set correct new angle	
	aa= thold + dtcorr;

end
function [THICK,XTHICK, CAMBR,XCAMBR] = tccalc(X,XP,Y,YP,S,N)
% ////--------------------------------------
% //     tccalc function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------


% C---------------------------------------------------------------
% C     Calculates max thickness and camber at airfoil pos
% C
% C     Note: this routine does not find the maximum camber or 
% C           thickness exactly as it only looks at discrete pos
% C
% C     Input:
% C       N      number of pos
% C       X(.)   shape coordinate po arrays
% C       Y(.)
% C
% C     Output:
% C       THICK  max thickness
% C       CAMBR  max camber
% C---------------------------------------------------------------
      [SLE]= lefind(X,XP,Y,YP,S,N);
      XLE = seval(SLE,X,XP,S,N);
      YLE = seval(SLE,Y,YP,S,N);
      XTE = 0.5*(X(1)+X(N));
      YTE = 0.5*(Y(1)+Y(N));
      CHORD = sqrt ((XTE-XLE)^2 + (YTE-YLE)^2);

% C---- set unit chord-line vector
      DXC = (XTE-XLE) / CHORD;
      DYC = (YTE-YLE) / CHORD;

      THICK = 0.;
      XTHICK = 0.;
      CAMBR = 0.;
      XCAMBR = 0.;
% C
% C---- go over each po, finding the y-thickness and camber
      for I=1: N
        XBAR = (X(I)-XLE)*DXC + (Y(I)-YLE)*DYC;
        YBAR = (Y(I)-YLE)*DXC - (X(I)-XLE)*DYC;
% C
% C------ set po on the opposite side with the same chord x value
        [SOPP]= sopps(S(I), X,XP,Y,YP,S,N, SLE);
        XOPP = seval(SOPP,X,XP,S,N);
        YOPP = seval(SOPP,Y,YP,S,N);

        YBAROP = (YOPP-YLE)*DXC - (XOPP-XLE)*DYC;

        YC = 0.5*(YBAR+YBAROP);
        YT =  abs(YBAR-YBAROP);

        if (abs(YC) > abs(CAMBR)) 
         CAMBR = YC;
         XCAMBR = XOPP;
        end
        if (abs(YT) > abs(THICK)) 
         THICK = YT;
         XTHICK = XOPP;
        end
      end

end
function sopp=sopps(  si,  x,  xp,  y,  yp,  s, n,  sle)
% ////--------------------------------------
% //     sopps function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% 	//      dimension x(*),xp(*),y(*),yp(*),s(*)
% 	//--------------------------------------------------
% 	//     calculates arc length sopp of point 
% 	//     which is opposite of point si, on the 
% 	//     other side of the airfoil baseline
% 	//--------------------------------------------------
% 	 chord, slen, dxc, dyc, sfrac;
% 	 xi, yi, xbar, xopp, yopp, xoppd, yoppd;
% 	 res, resd, dsopp;
% 	int in, inopp;
% 	//---- reference length for testing convergence
	slen = s(n) - s(1);
	
% 	//---this fails miserably with sharp le foils, tsk,tsk,tsk hhy 4/24/01
% 	//---- set baseline vector normal to surface at le point
% 	//      dxc = -deval(sle,y,yp,s,n)
% 	//      dyc =  deval(sle,x,xp,s,n)
% 	//      dsc = sqrt(dxc**2 + dyc**2)
% 	//      dxc = dxc/dsc
% 	//      dyc = dyc/dsc
% 	
% 	//---rational alternative 4/24/01 hhy
	xle = seval(sle,x,xp,s,n);
	yle = seval(sle,y,yp,s,n);
	xte = 0.5*(x(1)+x(n));
	yte = 0.5*(y(1)+y(n));
	chord = sqrt((xte-xle)*(xte-xle) + (yte-yle)*(yte-yle));
% 	//---- set unit chord-line vector
	dxc = (xte-xle) / chord;
	dyc = (yte-yle) / chord;
	
	if(si<sle) 
		in = 1;
		inopp = n;
	
    else
		in = n;
		inopp = 1;
    end
	sfrac = (si-sle)/(s(in)-sle);
	sopp = sle + sfrac*(s(inopp)-sle);
	
	if(abs(sfrac) <= 1.0e-5) 
		sopp = sle;
		return;
    end
	
% 	//---- xbar = x coordinate in chord-line axes
	xi  = seval(si , x,xp,s,n);
	yi  = seval(si , y,yp,s,n);
	xle = seval(sle, x,xp,s,n);
	yle = seval(sle, y,yp,s,n);
	xbar = (xi-xle)*dxc + (yi-yle)*dyc;
	
% 	//---- converge on exact opposite point with same xbar value
	bFound = false;
	for (iter=1:12)
		xopp  = seval(sopp,x,xp,s,n);
		yopp  = seval(sopp,y,yp,s,n);
		xoppd = deval(sopp,x,xp,s,n);
		yoppd = deval(sopp,y,yp,s,n);
		
		res  = (xopp -xle)*dxc + (yopp -yle)*dyc - xbar;
		resd =  xoppd     *dxc +  yoppd     *dyc;
		
		if(abs(res)/slen < 1.0e-5) 
			bFound = true;
			break%;//go to 305
        end
		if(resd == 0.0) 
			bFound = false;
			break%;//go to 303
        end
		dsopp = -res/resd;
		sopp =sopp+dsopp;
		
		if(abs(dsopp)/slen < 1.0e-5) 
			bFound = true;
			break%;//go to 305
        end
    end
	if(bFound==false)
% 		//303		write(*,*) 'sopps: opposite-point location failed. continuing...'
		sopp = sle + sfrac*(s(inopp)-sle);
    end
% 	// 305  continue
	
end
function p=abcopy(p)
% ////--------------------------------------
% //     abcopy function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% 	int i;
% 	if(p.nb<=1) 
% 		CString str;
% 		str.Format("abcopy: buffer airfoil not available");
% 		AfxMessageBox(str, MB_OK);
		
% 		return false;
% 	}
% 	else if(nb>IQX-2) {
% 		CString str1, str2;
% 		str1.Format("Maximum number of panel nodes  : %d\r\p.n",IQX-2);
% 		str2.Format("Number of buffer airfoil points: %d\r\p.n",nb);
% 		str2+="Current airfoil cannot be set\r\p.n";
% 		str2+="Try executing PANE at top level instead";
% 		str1+=str2;
% 		AfxMessageBox(str1, MB_OK);
% 		return false;
% 	}
	if(p.n~=p.nb) 
        p.lblini = false;
    end
	p.n = p.nb;
	for (i=1:p.n)
		p.x(i) = p.xb(i);
		p.y(i) = p.yb(i);
    end
	p.lgsame = true;
	
	if(p.lbflap) 
		p.xof = p.xbf;
		p.yof = p.ybf;
		p.lflap = true;
    end
	
% 	//---- strip out doubled points
	i = 1;
	
	while (i<p.n)
		i=i+1;
		if(p.x(i-1)==p.x(i) && p.y(i-1)==p.y(i)) 
			for (j=i:p.n-1)
				p.x(j) = p.x(j+1);
				p.y(j) = p.y(j+1);
            end
			p.n = p.n-1;
        end
    end

	p.s=scalc(p.x,p.y,p.n);
	p.xp=segspl(p.x,p.s,p.n);
	p.yp=segspl(p.y,p.s,p.n);
	[p.nx,p.ny]=ncalc(p.x,p.y,p.s,p.n);
	p.sle=lefind(p.x,p.xp,p.y,p.yp,p.s,p.n);
	p.xle = seval(p.sle,p.x,p.xp,p.s,p.n);
	p.yle = seval(p.sle,p.y,p.yp,p.s,p.n);
	p.xte = 0.5*(p.x(1)+p.x(p.n));
	p.yte = 0.5*(p.y(1)+p.y(p.n));
	p.chord  = sqrt( (p.xte-p.xle)*(p.xte-p.xle) + (p.yte-p.yle)*(p.yte-p.yle) );
	p=tecalc(p);
	p=apcalc(p);
	
	p.lgamu = false;
	p.lqinu = false;
	p.lwake = false;
	p.lqaij = false;
	p.ladij = false;
	p.lwdij = false;
	p.lipan = false;
	p.lvconv = false;
% //	lscini = false;
	
% 	//   write(*,1200) p.n
% 	// 1200 format(/' current airfoil nodes set from buffer airfoil nodes (', i4,' )')
	
% 	return true;
% }
end
function [xn,yn]=ncalc( x,  y,  s, n)
% ////--------------------------------------
% //     ncalc function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% {
% 	 sx, sy, smod;
% 	int i;	
	if(n<=1)
        return;
    end
	xn=segspl(x,s,n);
	yn=segspl(y,s,n);
	for (i=1:n)
		sx =  yn(i);
		sy = -xn(i);
		smod = sqrt(sx*sx + sy*sy);
		xn(i) = sx/smod;
		yn(i) = sy/smod;
    end
	
% 	//---- average normal vectors at corner points
	for (i=1:n-1)
		if(s(i) == s(i+1)) 
			sx = 0.5*(xn(i) + xn(i+1));
			sy = 0.5*(yn(i) + yn(i+1));
			smod = sqrt(sx*sx + sy*sy);
			xn(i)   = sx/smod;
			yn(i)   = sy/smod;
			xn(i+1) = sx/smod;
			yn(i+1) = sy/smod;
        end
    end
	
% 	return true;
% }
end
function p=tecalc(p)
% ////--------------------------------------
% //     tecalc function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% //-------------------------------------------
% //     calculates total and projected te 
% //     areas and te panel strengths.
% 	//-------------------------------------------
% 	
% 	 p.scs, p.sds;
% 	//---- set te base vector and te bisector components
	 p.dxte = p.x(1) - p.x(p.n);
	 p.dyte = p.y(1) - p.y(p.n);
	 p.dxs = 0.5*(-p.xp(1) + p.xp(p.n));
	 p.dys = 0.5*(-p.yp(1) + p.yp(p.n));
	
% 	//---- normal and streamwise projected te gap areas
	p.ante = p.dxs*p.dyte - p.dys*p.dxte;
	p.aste = p.dxs*p.dxte + p.dys*p.dyte;
	
% 	//---- total te gap area
	p.dste = sqrt(p.dxte*p.dxte + p.dyte*p.dyte);
	
	p.sharp = p.dste < 0.0001*p.chord;
	
	if(p.sharp) 
		p.scs = 1.0;
		p.sds = 0.0;
	
    else
		p.scs = p.ante/p.dste;
		p.sds = p.aste/p.dste;
    end
	
% 	//---- te panel source and vorticity strengths
	p.sigte = 0.5*(p.gam(1) - p.gam(p.n))*p.scs;
	p.gamte = -.5*(p.gam(1) - p.gam(p.n))*p.sds;
	
% //	sigte_a = 0.5*(gam_a(1) - gam_a(p.n))*p.scs;
% //	gamte_a = -.5*(gam_a(1) - gam_a(p.n))*p.sds;
	
% 	return true;
	
end
function p=apcalc(p)
% ////--------------------------------------
% //     apcalc function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% 	double p.sx, p.sy;
% 	int i, ip;
% 	//---- set angles of airfoil panels
	for (i=1:p.n-1)
		sx = p.x(i+1) - p.x(i);
		sy = p.y(i+1) - p.y(i);
		if(sx==0.0 && sy==0.0) 
            p.apanel(i) = atan2(-p.ny(i), -p.nx(i));
        else
            p.apanel(i) = atan2(sx, -sy );
        end
    end
		
	
	
% 	//---- te panel
	i = p.n;
	ip = 1;
	if(p.sharp) 
        p.apanel(i) = pi;
    else
		sx = p.x(ip) - p.x(i);
		sy = p.y(ip) - p.y(i);
		p.apanel(i) = atan2( -sx , sy ) + pi;
    end
	
% 	return true;
end
function [imax,amax]=cang( x,  y,  n)
% ////--------------------------------------
% //     cang function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% {
% //-------------------------------------------------------------------
% 	 dx1, dx2, dy1, dy2, crossp, angl;	
	amax = 0.0;
	imax = 1;
	
% 	//---- go over each point, calculating corner angle
	for ( i=2:n-1)
		dx1 = x(i) - x(i-1);
		dy1 = y(i) - y(i-1);
		dx2 = x(i) - x(i+1);
		dy2 = y(i) - y(i+1);
		
% 		//------ allow for doubled points
		if(dx1==0.0 && dy1==0.0) 
			dx1 = x(i) - x(i-2);
			dy1 = y(i) - y(i-2);
        end
		if(dx2==0.0 && dy2==0.0) 
			dx2 = x(i) - x(i+2);
			dy2 = y(i) - y(i+2);
        end
		
		crossp = (dx2*dy1 - dy2*dx1)...
			/ sqrt((dx1*dx1 + dy1*dy1) * (dx2*dx2 + dy2*dy2));
		angl = asin(crossp)*(180.0/pi);
		if(abs(angl) > abs(amax)) 
			amax = angl;
			imax = i;
        end
    end
% 	return true;
end
function p=specal(p)
% ////--------------------------------------
% //     specal function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% {
% //-----------------------------------
% //     converges to specified alpha.
% //-----------------------------------
% 	double minf_clm, p.msq_clm, reinf_clm;
% 	double clm, dclm, clm1;
% 	int i, irlx, itcl;
% 
% 	//---- calculate surface vorticity distributions for alpha = 0, 90 degrees
	if(p.lgamu==false || p.lqaij==false)
        p=ggcalc(p);
    end
	
	p.cosa = cos(p.alfa);
	p.sina = sin(p.alfa);
	
% 	//---- superimpose suitably weighted  alpha = 0, 90  distributions
	for (i=1:p.n)
		p.gam(i)   =  p.cosa*p.gamu(i,1) + p.sina*p.gamu(i,2);
		p.gam_a(i) = -p.sina*p.gamu(i,1) + p.cosa*p.gamu(i,2);
    end
	p.psio = p.cosa*p.gamu(p.n+1,1) + p.sina*p.gamu(p.n+1,2);
	
	p=tecalc(p);
	p=qiset(p);
	
% 	//---- set initial guess for the newton variable clm
	clm = 1.0;

% 	//---- set corresponding  m(clm), re(clm)
	[p,minf_clm,reinf_clm]=mrcl(clm,p);
	p=comset(p);
	
% 	//---- set corresponding cl(m)
	p=clcalc(p.xcmref,p.ycmref,p);
% 	//---- iterate on clm
	p.bConv = false;
	for (itcl=1:20)
		
		p.msq_clm = 2.0*p.minf*minf_clm;
		dclm = (p.cl - clm)/(1.0 - p.cl_msq*p.msq_clm);
		
		clm1 = clm;
		rlx = 1.0;
		
% 		//------ under-relaxation loop to avoid driving m(cl) above 1
		for (irlx=1:12)
			
			clm = clm1 + rlx*dclm;
			
% 			//-------- set new freestream mach m(clm)
			[p,minf_clm,reinf_clm]=mrcl(clm,p);
			
% 			//-------- if mach is ok, go do next newton iteration
			if(p.matyp==1 || p.minf==0.0 || minf_clm~=0.0) 
                break
            end
			rlx = 0.5*rlx;
        end

		
% 		//------ set new cl(m)
		p=comset(p);
		p=clcalc(p.xcmref,p.ycmref,p);
		
		if(abs(dclm)<=1.0e-6) 
			p.bConv = true;
			break;
        end
    end
	if(p.bConv==false)		
% //		AfxMessageBox("Specal:  MInf convergence failed", MB_OK);
% 		CString str;
		disp('Specal:  MInf convergence failed');
% 		if(m_bTrace)pXFile->WriteString(str);
		return;
%         end
    end
	
% 	//---- set final mach, cl, cp distributions, and hinge moment
	[p,minf_cl,reinf_cl]=mrcl(p.cl,p);
	p=comset(p);
	p=clcalc(p.xcmref,p.ycmref,p);
% /*
% 	if (!cpcalc(n,qinv,qinf,p.minf,cpi)) return false;// no need to carry on
% 	if(lvisc) {
% 		if(!cpcalc(n+nw,qvis,qinf,p.minf,cpv)) return false;// no need to carry on
%         end
% 		if(!cpcalc(n+nw,qinv,qinf,p.minf,cpi)) return false;// no need to carry on
% 	}
% 	else   if (!cpcalc(n,qinv,qinf,p.minf,cpi)) return false;// no need to carry on
% */
% 	p.cpi=cpcalc(p.n,p.qinv,p.qinf,p.minf);
%     p.cpi
% 	if(p.lvisc) 
% 	{
% 		cpcalc(n+nw,qvis,qinf,p.minf,cpv);
% 		cpcalc(n+nw,qinv,qinf,p.minf,cpi);
% 	}
% 	else 
if (p.lvisc==false)
 		p.cpi=cpcalc(p.n,p.qinv,p.qinf,p.minf);
end
	
%     end
	if(p.lflap)
        p=mhinge(p);
    end
% 	//Added arcds to get inviscid q after viscous calculation
	for (i=1:p.n)
		p.qgamm(i) = p.gam(i);
    end
% 	// end arcds addition

% 	return true;
end
function p=ggcalc(p)
% ////--------------------------------------
% //     ggcalc function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------
% {
% //--------------------------------------------------------------
% //     calculates two surface vorticity (gamma) distributions
% //     for alpha = 0, 90  degrees.  these are superimposed
% //     in specal or speccl for specified alpha or cl.
% //--------------------------------------------------------------
	% int i,j, iu;
	% double psi, psi_n, p.psiinf, res, res1, res2, p.ag1, p.ag2;
	% double p.abis, p.cbis, p.sbis, p.ds1, p.ds2, p.dsmin;
	% double p.xbis, p.ybis, p.qbis;
	% p.cosa = cos(alfa);
	% sina = sin(alfa);
% 
% //---- distance of internal control point ahead of sharp te
% //-    (fraction of smaller panel length adjacent to te)
	bwt = 0.1;
	
	disp('calculating unit vorticity distributions ...');
% 	CString str;
% 	str.Format(" Calculating unit vorticity distributions ...\r\p.n");
% 	if(m_bTrace)pXFile->WriteString(str);

	for(i=1:p.n)
		p.gam(i) = 0.0;
		p.gamu(i,1) = 0.0;
		p.gamu(i,2) = 0.0;
    end
	p.psio = 0.0;
	
% 	//---- set up matrix system for  psi = p.psio  on airfoil surface.
% 	//-    the unknowns are (dgamma)i and dpsio.
	for (i=1:p.n)
		
% 		//------ calculate psi and dpsi/dgamma array for current node
		[psi,psi_n,p]=psilin(i,p.x(i),p.y(i),p.nx(i),p.ny(i),false,true,p);
		
		p.psiinf = p.qinf*(p.cosa*p.y(i) - p.sina*p.x(i));
		
% 		//------ res1 = psi( 0) - p.psio
% 		//------ res2 = psi(90) - p.psio
		res1 =  p.qinf*p.y(i);
		res2 = -p.qinf*p.x(i);
		
% 		//------ dres/dgamma
		for (j=1:p.n)
			p.aij(i,j) = p.dzdg(j);
        end
		
		for (j=1:p.n)
			p.bij(i,j) = -p.dzdm(j);
        end
		
% 		//------ dres/dpsio
		p.aij(i,p.n+1) = -1.0;
		
		p.gamu(i,1) = -res1;
		p.gamu(i,2) = -res2;
    end
	
% 	//---- set kutta condition
% 	//-    res = p.gam(1) + p.gam(p.n)
	res = 0.0;
	
	for (j=1:p.n+1) 
        p.aij(p.n+1,j) = 0.0;
    end
	
	
	p.aij(p.n+1,1) = 1.0;
	p.aij(p.n+1,p.n) = 1.0;
	
	p.gamu(p.n+1,1) = -res;
	p.gamu(p.n+1,2) = -res;
	
% 	//---- set up kutta condition (no direct source influence)
	for (j=1:p.n)
        p.bij(p.n+1,j) = 0.0;
    end
	
	if(p.sharp)
% 		//----- set zero internal velocity in te corner 
		
% 		//----- set te bisector angle
		p.ag1 = atan2(-p.yp(1),-p.xp(1)    );
		p.ag2 = atanc( p.yp(p.n), p.xp(p.n),p.ag1);
		p.abis = 0.5*(p.ag1+p.ag2);
		p.cbis = cos(p.abis);
		p.sbis = sin(p.abis);
		
% 		//----- minimum panel length adjacent to te
		p.ds1 = sqrt((p.x(1)-p.x(2)  )*(p.x(1)-p.x(2)  ) + (p.y(1)-p.y(2)  )*(p.y(1)-p.y(2)  ));
		p.ds2 = sqrt((p.x(p.n)-p.x(p.n-1))*(p.x(p.n)-p.x(p.n-1)) + (p.y(p.n)-p.y(p.n-1))*(p.y(p.n)-p.y(p.n-1)));
		p.dsmin = min( p.ds1 , p.ds2 );
		
% 		//----- control point on bisector just ahead of te point
		p.xbis = p.xte - bwt*p.dsmin*p.cbis;
		p.ybis = p.yte - bwt*p.dsmin*p.sbis;
		
% 		//----- set velocity component along bisector line
		[psi,p.qbis,p]=psilin(0,p.xbis,p.ybis,-p.sbis,p.cbis,false,true,p);
		
		res = p.qbis;
		
% 		//----- dres/dgamma
		for (j=1:p.n) 
            p.aij(p.n,j) = dqdg(j);
        end
		
		
% 		//----- -dres/dmass
		for (j=1:p.n) 
			p.bij(p.n,j) = -dqdm(j);
        end
		
% 		//----- dres/dpsio
		p.aij(p.n,p.n+1) = 0.0;
		
% 		//----- -dres/duinf
		p.gamu(p.n,1) = -p.cbis;
		
% 		//----- -dres/dvinf
		p.gamu(p.n,2) = -p.sbis;
		
    end
	
% 	//---- lu-factor coefficient matrix p.aij
	[p.aij,p.aijpiv]=ludcmp(p.n+1,p.aij);
	p.lqaij = true;
	
% 	//---- solve system for the two vorticity distributions
% 	double bbb(IQX);
%     IQX=360;
% 	for (iu=0:IQX) bbb(iu) = p.gamu(iu,1);//arcds : create a dummy array
      p.gamu(:,1)= baksub(p.n+1,p.aij,p.aijpiv,p.gamu(:,1));
      p.gamu(:,2)= baksub(p.n+1,p.aij,p.aijpiv,p.gamu(:,2));
% 	baksub(p.n+1, p.aij, aijpiv, bbb);
% 	for (iu=0; iu<IQX; iu++) p.gamu(iu,1) = bbb(iu) ;
	
% 	for (iu=0; iu<IQX; iu++) bbb(iu) = p.gamu(iu,2);//arcds : create a dummy array
% 	baksub(p.n+1, p.aij, aijpiv, bbb);
% 	for (iu=0; iu<IQX; iu++) p.gamu(iu,2) = bbb(iu) ;

	
% 	//---- set inviscid alpha=0,90 surface speeds for this geometry
	for (i=1:p.n)
		p.qinvu(i,1) = p.gamu(i,1);
		p.qinvu(i,2) = p.gamu(i,2);
    end
	
	p.lgamu = true;
	
% 	return true;

 
end
function [psi,psi_ni,p]=psilin(i,xi,yi,nxi,nyi,geolin,siglin,p)
% ////--------------------------------------
% //     psilin function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% function [PSI,PSI_NI,DZDG,DZDN,DQDG,DZDM,DQDM,Z_QINF,Z_ALFA,Z_QDOF0,Z_QDOF1,Z_QDOF2,Z_QDOF3,QTAN1,QTAN2,QTANM,SX,SY]=PSILIN(I,XI,YI,NXI,NYI,GEOLIN,SIGLIN,S,N,SHARP,ANTE,ASTE,DSTE,X,Y,APANEL,NX,NY,SIG,GAMU,GAM,ALFA)
% 
% % C-----------------------------------------------------------------------
% % C     Calculates current streamfunction Psi at panel node or wake node
% % C     I due to freestream and all bound vorticity Gam on the airfoil. 
% % C     Sensitivities of Psi with respect to alpha (Z_ALFA) and inverse
% % C     Qspec DOFs (Z_QDOF0,Z_QDOF1) which influence Gam in inverse cases.
% % C     Also calculates the sensitivity vector dPsi/dGam (DZDG).
% % C
% % C     If SIGLIN=True, then Psi includes the effects of the viscous
% % C     source distribution Sig and the sensitivity vector dPsi/dSig
% % C     (DZDM) is calculated.
% % C
% % C     If GEOLIN=True, then the geometric sensitivity vector dPsi/dn
% % C     is calculated, where p.n is the normal motion of the jth node.
% % C
% % C          Airfoil:  1   < I < N
% % C          Wake:     N+1 < I < N+NW
% % C-----------------------------------------------------------------------
% % C
% % C---- distance tolerance for determining if two pos are the same
% COSA=cos(ALFA);
% SINA=sin(ALFA);
% 	QOPI = 0.25/pi;
%     HOPI = 0.50/pi;
%  QINF = 1.0;
%  LIMAGE=false;
%        YIMAGE = -10.0;
% 
% 
%      SEPS = (S(N)-S(1)) * 0.00001;
%       IO = I;
% 
%   %    COSA = COS(ALFA);
%   %    SINA = SIN(ALFA);
% 
%       for JO=1: N
%         DZDG(JO) = 0.0;
%         DZDN(JO) = 0.0;
%         DQDG(JO) = 0.0;
%       end
% 
%       for JO=1: N
%         DZDM(JO) = 0.0;
%         DQDM(JO) = 0.0;
%       end
% 
%       Z_QINF = 0.;
%       Z_ALFA = 0.;
%       Z_QDOF0 = 0.;
%       Z_QDOF1 = 0.;
%       Z_QDOF2 = 0.;
%       Z_QDOF3 = 0.;
% 
%       PSI    = 0.;
%       PSI_NI = 0.;
% 
%       QTAN1 = 0.;
%       QTAN2 = 0.;
%       QTANM = 0.;
% 
%       if (SHARP~=0)
%        SCS = 1.0;
%        SDS = 0.0;
%       else
%        SCS = ANTE/DSTE;
%        SDS = ASTE/DSTE;
%       end
% 
%       %DO 10 
%       for JO=1: N
%         JP = JO+1;
% 
%         JM = JO-1;
%         JQ = JP+1;
% 
%         if (JO==1)
%          JM = JO;
%         elseif (JO==N-1) 
%          JQ = JP;
%         elseif(JO==N) 
%          JP = 1;
%          if((X(JO)-X(JP))^2 + (Y(JO)-Y(JP))^2 < SEPS^2) %GO TO 12
%              
%              
%                    PSI = PSI + QINF*(COSA*YI - SINA*XI);
% % C
% % C---- dPsi/dn
%       PSI_NI = PSI_NI + QINF*(COSA*NYI - SINA*NXI);
% 
%       QTAN1 = QTAN1 + QINF*NYI;
%       QTAN2 = QTAN2 - QINF*NXI;
% % C
% % C---- dPsi/dQinf
%       Z_QINF = Z_QINF + (COSA*YI - SINA*XI);
% % C
% % C---- dPsi/dalfa
%       Z_ALFA = Z_ALFA - QINF*(SINA*YI + COSA*XI);
% % C
%       if (LIMAGE==false) 
%       return
%       end
%       
%       for JO=1: N
%         JP = JO+1;
% 
%         JM = JO-1;
%         JQ = JP+1;
% 
%         if(JO==1)
%          JM = JO;
%         elseif (JO==N-1) 
%          JQ = JP;
%         elseif (JO==N) 
%          JP = 1;
%          if((X(JO)-X(JP))^2 + (Y(JO)-Y(JP))^2 < SEPS^2) 
%            return
%          end
%         end
% 
%         DSO = sqrt((X(JO)-X(JP))^2 + (Y(JO)-Y(JP))^2);
% % C
% % C------ skip null panel
%         if (DSO == 0.0) 
%             continue
%         end
% 
%         DSIO = 1.0 / DSO;
% 
%         APAN = pi - APANEL(JO) + 2.0*ALFA;
% 
%         XJO = X(JO) + 2.0*(YIMAGE+Y(JO))*SINA;
%         YJO = Y(JO) - 2.0*(YIMAGE+Y(JO))*COSA;
%         XJP = X(JP) + 2.0*(YIMAGE+Y(JP))*SINA;
%         YJP = Y(JP) - 2.0*(YIMAGE+Y(JP))*COSA;
% 
%         RX1 = XI - XJO;
%         RY1 = YI - YJO;
%         RX2 = XI - XJP;
%         RY2 = YI - YJP;
% 
%         SX = (XJP - XJO) * DSIO;
%         SY = (YJP - YJO) * DSIO;
% 
%         X1 = SX*RX1 + SY*RY1;
%         X2 = SX*RX2 + SY*RY2;
%         YY = SX*RY1 - SY*RX1;
% 
%         RS1 = RX1*RX1 + RY1*RY1;
%         RS2 = RX2*RX2 + RY2*RY2;
% % C
% % C------ set reflection flag SGN to avoid branch problems with arctan
%         if(IO>=1 && IO~=N) 
% %C------- no problem on airfoil surface
%          SGN = 1.0;
%         else
% % C------- make sure arctan falls between  -/+  Pi/2
%          SGN = 1*sign(YY);
%         end
% 
% % C------ set log(r^2) and arctan(p.x/p.y), correcting for reflection if any
%         G1 = log(RS1);
%         T1 = atan2(SGN*X1,SGN*YY) + (0.5 - 0.5*SGN)*pi;
% 
%         G2 = log(RS2);
%         T2 = atan2(SGN*X2,SGN*YY) + (0.5 - 0.5*SGN)*pi;
% 
%         X1I = SX*NXI + SY*NYI;
%         X2I = SX*NXI + SY*NYI;
%         YYI = SX*NYI - SY*NXI;
% 
%         if (GEOLIN==true)
%          NXO = NX(JO);
%          NYO = NY(JO);
%          NXP = NX(JP);
%          NYP = NY(JP);
% 
%          X1O =-((RX1-X1*SX)*NXO + (RY1-X1*SY)*NYO)*DSIO-(SX*NXO+SY*NYO);
%          X1P = ((RX1-X1*SX)*NXP + (RY1-X1*SY)*NYP)*DSIO;
%          X2O =-((RX2-X2*SX)*NXO + (RY2-X2*SY)*NYO)*DSIO;
%          X2P = ((RX2-X2*SX)*NXP + (RY2-X2*SY)*NYP)*DSIO-(SX*NXP+SY*NYP);
%          YYO = ((RX1+X1*SY)*NYO - (RY1-X1*SX)*NXO)*DSIO-(SX*NYO-SY*NXO);
%          YYP =-((RX1-X1*SY)*NYP - (RY1+X1*SX)*NXP)*DSIO;
%         end
% 
%         if (JO==N) 
%             break
%         end
% 
%         if (SIGLIN==true) 
% % C
% % C------- set up midpoint quantities
%          X0 = 0.5*(X1+X2);
%          RS0 = X0*X0 + YY*YY;
%          G0 = log(RS0);
%          T0 = atan2(SGN*X0,SGN*YY) + (0.5 - 0.5*SGN)*pi;
% % C
% % C------- calculate source contribution to Psi  for  1-0  half-panel
%          DXINV = 1.0/(X1-X0);
%          PSUM = X0*(T0-APAN) - X1*(T1-APAN) + 0.5*YY*(G1-G0);
%          PDIF = ((X1+X0)*PSUM + RS1*(T1-APAN) - RS0*(T0-APAN)    + (X0-X1)*YY) * DXINV;
% 
%          PSX1 =  -(T1-APAN);
%          PSX0 =    T0-APAN;
%          PSYY =  0.5*(G1-G0);
% 
%          PDX1 = ((X1+X0)*PSX1 + PSUM + 2.0*X1*(T1-APAN) - PDIF) * DXINV;
%          PDX0 = ((X1+X0)*PSX0 + PSUM - 2.0*X0*(T0-APAN) + PDIF) * DXINV;
%          PDYY = ((X1+X0)*PSYY + 2.0*(X0-X1 + YY*(T1-T0))      ) * DXINV;
% 
%          DSM = sqrt((X(JP)-X(JM))^2 + (Y(JP)-Y(JM))^2);
%          DSIM = 1.0/DSM;
% 
% % CCC      SIG0 = (SIG(JP) - SIG(JO))*DSIO
% % CCC      SIG1 = (SIG(JP) - SIG(JM))*DSIM
% % CCC      SSUM = SIG0 + SIG1
% % CCC      SDIF = SIG0 - SIG1
% % C
%          SSUM = (SIG(JP) - SIG(JO))*DSIO + (SIG(JP) - SIG(JM))*DSIM;
%          SDIF = (SIG(JP) - SIG(JO))*DSIO - (SIG(JP) - SIG(JM))*DSIM;
% 
%          PSI = PSI + QOPI*(PSUM*SSUM + PDIF*SDIF);
% % C
% % C------- dPsi/dm
%          DZDM(JM) = DZDM(JM) + QOPI*(-PSUM*DSIM + PDIF*DSIM);
%          DZDM(JO) = DZDM(JO) + QOPI*(-PSUM*DSIO - PDIF*DSIO);
%          DZDM(JP) = DZDM(JP) + QOPI*( PSUM*(DSIO+DSIM)     + PDIF*(DSIO-DSIM));
% % C
% % C------- dPsi/dni
%          PSNI = PSX1*X1I + PSX0*(X1I+X2I)*0.5 + PSYY*YYI;
%          PDNI = PDX1*X1I + PDX0*(X1I+X2I)*0.5 + PDYY*YYI;
%          PSI_NI = PSI_NI + QOPI*(PSNI*SSUM + PDNI*SDIF);
% 
%          QTANM = QTANM + QOPI*(PSNI*SSUM + PDNI*SDIF);
% 
%          DQDM(JM) = DQDM(JM) + QOPI*(-PSNI*DSIM + PDNI*DSIM);
%          DQDM(JO) = DQDM(JO) + QOPI*(-PSNI*DSIO - PDNI*DSIO);
%          DQDM(JP) = DQDM(JP) + QOPI*( PSNI*(DSIO+DSIM) + PDNI*(DSIO-DSIM));
% 
% 
% %C------- calculate source contribution to Psi  for  0-2  half-panel
%          DXINV = 1.0/(X0-X2);
%          PSUM = X2*(T2-APAN) - X0*(T0-APAN) + 0.5*YY*(G0-G2);
%          PDIF = ((X0+X2)*PSUM + RS0*(T0-APAN) - RS2*(T2-APAN) + (X2-X0)*YY) * DXINV;
% 
%          PSX0 =  -(T0-APAN);
%          PSX2 =    T2-APAN;
%          PSYY =  0.5*(G0-G2);
% 
%          PDX0 = ((X0+X2)*PSX0 + PSUM + 2.0*X0*(T0-APAN) - PDIF) * DXINV;
%          PDX2 = ((X0+X2)*PSX2 + PSUM - 2.0*X2*(T2-APAN) + PDIF) * DXINV;
%          PDYY = ((X0+X2)*PSYY + 2.0*(X2-X0 + YY*(T0-T2))      ) * DXINV;
% 
%          DSP = sqrt((X(JQ)-X(JO))^2 + (Y(JQ)-Y(JO))^2);
%          DSIP = 1.0/DSP;
% % 
% % CCC         SIG2 = (SIG(JQ) - SIG(JO))*DSIP
% % CCC         SIG0 = (SIG(JP) - SIG(JO))*DSIO
% % CCC         SSUM = SIG2 + SIG0
% % CCC         SDIF = SIG2 - SIG0
% % C
%          SSUM = (SIG(JQ) - SIG(JO))*DSIP + (SIG(JP) - SIG(JO))*DSIO;
%          SDIF = (SIG(JQ) - SIG(JO))*DSIP - (SIG(JP) - SIG(JO))*DSIO;
% 
%          PSI = PSI + QOPI*(PSUM*SSUM + PDIF*SDIF);
% 
% %C------- dPsi/dm
%          DZDM(JO) = DZDM(JO) + QOPI*(-PSUM*(DSIP+DSIO)      - PDIF*(DSIP-DSIO));
%          DZDM(JP) = DZDM(JP) + QOPI*( PSUM*DSIO - PDIF*DSIO);
%          DZDM(JQ) = DZDM(JQ) + QOPI*( PSUM*DSIP + PDIF*DSIP);
% 
% %C------- dPsi/dni
%          PSNI = PSX0*(X1I+X2I)*0.5 + PSX2*X2I + PSYY*YYI;
%          PDNI = PDX0*(X1I+X2I)*0.5 + PDX2*X2I + PDYY*YYI;
%          PSI_NI = PSI_NI + QOPI*(PSNI*SSUM + PDNI*SDIF);
% 
%          QTANM = QTANM + QOPI*(PSNI*SSUM + PDNI*SDIF);
% 
%          DQDM(JO) = DQDM(JO) + QOPI*(-PSNI*(DSIP+DSIO)          - PDNI*(DSIP-DSIO));
%          DQDM(JP) = DQDM(JP) + QOPI*( PSNI*DSIO - PDNI*DSIO);
%          DQDM(JQ) = DQDM(JQ) + QOPI*( PSNI*DSIP + PDNI*DSIP);
% 
%         end
% 
% %C------ calculate vortex panel contribution to Psi
%         DXINV = 1.0/(X1-X2);
%         PSIS = 0.5*X1*G1 - 0.5*X2*G2 + X2 - X1 + YY*(T1-T2);
%         PSID = ((X1+X2)*PSIS + 0.5*(RS2*G2-RS1*G1 + X1*X1-X2*X2))*DXINV;
% 
%         PSX1 = 0.5*G1;
%         PSX2 = -.5*G2;
%         PSYY = T1-T2;
% 
%         PDX1 = ((X1+X2)*PSX1 + PSIS - X1*G1 - PSID)*DXINV;
%         PDX2 = ((X1+X2)*PSX2 + PSIS + X2*G2 + PSID)*DXINV;
%         PDYY = ((X1+X2)*PSYY - YY*(G1-G2)         )*DXINV;
% 
%         GSUM1 = GAMU(JP,1) + GAMU(JO,1);
%         GSUM2 = GAMU(JP,2) + GAMU(JO,2);
%         GDIF1 = GAMU(JP,1) - GAMU(JO,1);
%         GDIF2 = GAMU(JP,2) - GAMU(JO,2);
% 
%         GSUM = GAM(JP) + GAM(JO);
%         GDIF = GAM(JP) - GAM(JO);
% 
%         PSI = PSI - QOPI*(PSIS*GSUM + PSID*GDIF);
% 
% %C------ dPsi/dGam
%         DZDG(JO) = DZDG(JO) - QOPI*(PSIS-PSID);
%         DZDG(JP) = DZDG(JP) - QOPI*(PSIS+PSID);
% 
% %C------ dPsi/dni
%         PSNI = PSX1*X1I + PSX2*X2I + PSYY*YYI;
%         PDNI = PDX1*X1I + PDX2*X2I + PDYY*YYI;
%         PSI_NI = PSI_NI - QOPI*(GSUM*PSNI + GDIF*PDNI);
% 
%         QTAN1 = QTAN1 - QOPI*(GSUM1*PSNI + GDIF1*PDNI);
%         QTAN2 = QTAN2 - QOPI*(GSUM2*PSNI + GDIF2*PDNI);
% 
%         DQDG(JO) = DQDG(JO) - QOPI*(PSNI - PDNI);
%         DQDG(JP) = DQDG(JP) - QOPI*(PSNI + PDNI);
% 
%         if (GEOLIN==true) 
% % C
% % C------- dPsi/dn
%          DZDN(JO) = DZDN(JO)- QOPI*GSUM*(PSX1*X1O + PSX2*X2O + PSYY*YYO) - QOPI*GDIF*(PDX1*X1O + PDX2*X2O + PDYY*YYO);
%          DZDN(JP) = DZDN(JP)- QOPI*GSUM*(PSX1*X1P + PSX2*X2P + PSYY*YYP) - QOPI*GDIF*(PDX1*X1P + PDX2*X2P + PDYY*YYP);
% %C------- dPsi/dP
%          Z_QDOF0 = Z_QDOF0         - QOPI*((PSIS-PSID)*QF0(JO) + (PSIS+PSID)*QF0(JP));
%          Z_QDOF1 = Z_QDOF1      - QOPI*((PSIS-PSID)*QF1(JO) + (PSIS+PSID)*QF1(JP));
%          Z_QDOF2 = Z_QDOF2       - QOPI*((PSIS-PSID)*QF2(JO) + (PSIS+PSID)*QF2(JP));
%          Z_QDOF3 = Z_QDOF3       - QOPI*((PSIS-PSID)*QF3(JO) + (PSIS+PSID)*QF3(JP));
%         end
% 
% 
%       end
%          
%       PSIG = 0.5*YY*(G1-G2) + X2*(T2-APAN) - X1*(T1-APAN);
%       PGAM = 0.5*X1*G1 - 0.5*X2*G2 + X2 - X1 + YY*(T1-T2);
% 
%       PSIGX1 = -(T1-APAN);
%       PSIGX2 =   T2-APAN;
%       PSIGYY = 0.5*(G1-G2);
%       PGAMX1 = 0.5*G1;
%       PGAMX2 = -.5*G2;
%       PGAMYY = T1-T2;
% 
%       PSIGNI = PSIGX1*X1I + PSIGX2*X2I + PSIGYY*YYI;
%       PGAMNI = PGAMX1*X1I + PGAMX2*X2I + PGAMYY*YYI;
% 
% %C---- TE panel source and vortex strengths
%       SIGTE1 = 0.5*SCS*(GAMU(JP,1) - GAMU(JO,1));
%       SIGTE2 = 0.5*SCS*(GAMU(JP,2) - GAMU(JO,2));
%       GAMTE1 = -.5*SDS*(GAMU(JP,1) - GAMU(JO,1));
%       GAMTE2 = -.5*SDS*(GAMU(JP,2) - GAMU(JO,2));
% 
%       SIGTE = 0.5*SCS*(GAM(JP) - GAM(JO));
%       GAMTE = -.5*SDS*(GAM(JP) - GAM(JO));
% 
% %C---- TE panel contribution to Psi
%       PSI = PSI + HOPI*(PSIG*SIGTE - PGAM*GAMTE);
% 
% %C---- dPsi/dGam
%       DZDG(JO) = DZDG(JO) - HOPI*PSIG*SCS*0.5;
%       DZDG(JP) = DZDG(JP) + HOPI*PSIG*SCS*0.5;
% 
%       DZDG(JO) = DZDG(JO) - HOPI*PGAM*SDS*0.5;
%       DZDG(JP) = DZDG(JP) + HOPI*PGAM*SDS*0.5;
% 
% %C---- dPsi/dni
%       PSI_NI = PSI_NI + HOPI*(PSIGNI*SIGTE - PGAMNI*GAMTE);
% 
%       QTAN1 = QTAN1 + HOPI*(PSIGNI*SIGTE1 - PGAMNI*GAMTE1);
%       QTAN2 = QTAN2 + HOPI*(PSIGNI*SIGTE2 - PGAMNI*GAMTE2);
% 
%       DQDG(JO) = DQDG(JO) - HOPI*(PSIGNI*0.5*SCS + PGAMNI*0.5*SDS);
%       DQDG(JP) = DQDG(JP) + HOPI*(PSIGNI*0.5*SCS + PGAMNI*0.5*SDS);
% 
%       if (GEOLIN==true) 
% 
% %C----- dPsi/dn
%        DZDN(JO) = DZDN(JO)  + HOPI*(PSIGX1*X1O + PSIGX2*X2O + PSIGYY*YYO)*SIGTE  - HOPI*(PGAMX1*X1O + PGAMX2*X2O + PGAMYY*YYO)*GAMTE;
%        DZDN(JP) = DZDN(JP) + HOPI*(PSIGX1*X1P + PSIGX2*X2P + PSIGYY*YYP)*SIGTE      - HOPI*(PGAMX1*X1P + PGAMX2*X2P + PGAMYY*YYP)*GAMTE;
% % C
% % C----- dPsi/dP
%        Z_QDOF0 = Z_QDOF0 + HOPI*PSIG*0.5*(QF0(JP)-QF0(JO))*SCS                + HOPI*PGAM*0.5*(QF0(JP)-QF0(JO))*SDS;
%        Z_QDOF1 = Z_QDOF1 + HOPI*PSIG*0.5*(QF1(JP)-QF1(JO))*SCS                + HOPI*PGAM*0.5*(QF1(JP)-QF1(JO))*SDS;
%        Z_QDOF2 = Z_QDOF2 + HOPI*PSIG*0.5*(QF2(JP)-QF2(JO))*SCS                + HOPI*PGAM*0.5*(QF2(JP)-QF2(JO))*SDS;
%        Z_QDOF3 = Z_QDOF3 + HOPI*PSIG*0.5*(QF3(JP)-QF3(JO))*SCS               + HOPI*PGAM*0.5*(QF3(JP)-QF3(JO))*SDS;
%       end
%          
% return
%          end
%         end
%              
%              
%              
%              
%         
% 
%         DSO = sqrt((X(JO)-X(JP))^2 + (Y(JO)-Y(JP))^2);
% 
% %C------ skip null panel
%         if (DSO == 0.0) %GO TO 10
%             continue
%         end
% 
%         DSIO = 1.0 / DSO;
% 
%         APAN = APANEL(JO);
% 
%         RX1 = XI - X(JO);
%         RY1 = YI - Y(JO);
%         RX2 = XI - X(JP);
%         RY2 = YI - Y(JP);
% 
%         SX = (X(JP) - X(JO)) * DSIO;
%         SY = (Y(JP) - Y(JO)) * DSIO;
% 
%         X1 = SX*RX1 + SY*RY1;
%         X2 = SX*RX2 + SY*RY2;
%         YY = SX*RY1 - SY*RX1;
% 
%         RS1 = RX1*RX1 + RY1*RY1;
%         RS2 = RX2*RX2 + RY2*RY2;
% 
% %C------ set reflection flag SGN to avoid branch problems with arctan
%         if (IO>=1 && IO~=N) 
% %C------- no problem on airfoil surface
%          SGN = 1.0;
%         else
% %C------- make sure arctan falls between  -/+  Pi/2
%          SGN = 1*sign(YY);
%         end
% % C
% % C------ set log(r^2) and arctan(p.x/p.y), correcting for reflection if any
%         if(IO~=JO && RS1>0.0) 
%          G1 = log(RS1);
%          T1 = atan2(SGN*X1,SGN*YY) + (0.5 - 0.5*SGN)*pi;
%         else
%          G1 = 0.0;
%          T1 = 0.0;
%         end
% 
%         if (IO~=JP && RS2>0.0) 
%          G2 = log(RS2);
%          T2 = atan2(SGN*X2,SGN*YY) + (0.5 - 0.5*SGN)*pi;
%         else
%          G2 = 0.0;
%          T2 = 0.0;
%         end
%        
% 
%         X1I = SX*NXI + SY*NYI;
%         X2I = SX*NXI + SY*NYI;
%         YYI = SX*NYI - SY*NXI;
% 
%         if (GEOLIN==true)
%          NXO = NX(JO);
%          NYO = NY(JO);
%          NXP = NX(JP);
%          NYP = NY(JP);
% 
%          X1O =-((RX1-X1*SX)*NXO + (RY1-X1*SY)*NYO)*DSIO-(SX*NXO+SY*NYO);
%          X1P = ((RX1-X1*SX)*NXP + (RY1-X1*SY)*NYP)*DSIO;
%          X2O =-((RX2-X2*SX)*NXO + (RY2-X2*SY)*NYO)*DSIO;
%          X2P = ((RX2-X2*SX)*NXP + (RY2-X2*SY)*NYP)*DSIO-(SX*NXP+SY*NYP);
%          YYO = ((RX1+X1*SY)*NYO - (RY1-X1*SX)*NXO)*DSIO-(SX*NYO-SY*NXO);
%          YYP =-((RX1-X1*SY)*NYP - (RY1+X1*SX)*NXP)*DSIO;
%         end
% 
%         if(JO==N) break
%         end
%         if (SIGLIN==true) 
% % C
% % C------- set up midpoint quantities
%          X0 = 0.5*(X1+X2);
%          RS0 = X0*X0 + YY*YY;
%          G0 = log(RS0);
%          T0 = atan2(SGN*X0,SGN*YY) + (0.5 - 0.5*SGN)*pi;
% 
% %C------- calculate source contribution to Psi  for  1-0  half-panel
%          DXINV = 1.0/(X1-X0);
%          PSUM = X0*(T0-APAN) - X1*(T1-APAN) + 0.5*YY*(G1-G0);
%          PDIF = ((X1+X0)*PSUM + RS1*(T1-APAN) - RS0*(T0-APAN)   + (X0-X1)*YY) * DXINV;
% 
%          PSX1 =  -(T1-APAN);
%          PSX0 =    T0-APAN;
%          PSYY =  0.5*(G1-G0);
% 
%          PDX1 = ((X1+X0)*PSX1 + PSUM + 2.0*X1*(T1-APAN) - PDIF) * DXINV;
%          PDX0 = ((X1+X0)*PSX0 + PSUM - 2.0*X0*(T0-APAN) + PDIF) * DXINV;
%          PDYY = ((X1+X0)*PSYY + 2.0*(X0-X1 + YY*(T1-T0))      ) * DXINV;
% 
%          DSM = sqrt((X(JP)-X(JM))^2 + (Y(JP)-Y(JM))^2);
%          DSIM = 1.0/DSM;
% % CCC      SIG0 = (SIG(JP) - SIG(JO))*DSIO
% % CCC      SIG1 = (SIG(JP) - SIG(JM))*DSIM
% % CCC      SSUM = SIG0 + SIG1
% % CCC      SDIF = SIG0 - SIG1
% % C
%          SSUM = (SIG(JP) - SIG(JO))*DSIO + (SIG(JP) - SIG(JM))*DSIM;
%          SDIF = (SIG(JP) - SIG(JO))*DSIO - (SIG(JP) - SIG(JM))*DSIM;
% 
%          PSI = PSI + QOPI*(PSUM*SSUM + PDIF*SDIF);
% 
% %C------- dPsi/dm
%          DZDM(JM) = DZDM(JM) + QOPI*(-PSUM*DSIM + PDIF*DSIM);
%          DZDM(JO) = DZDM(JO) + QOPI*(-PSUM*DSIO - PDIF*DSIO);
%          DZDM(JP) = DZDM(JP) + QOPI*( PSUM*(DSIO+DSIM)              + PDIF*(DSIO-DSIM));
% % C
% % C------- dPsi/dni
%          PSNI = PSX1*X1I + PSX0*(X1I+X2I)*0.5 + PSYY*YYI;
%          PDNI = PDX1*X1I + PDX0*(X1I+X2I)*0.5 + PDYY*YYI;
%          PSI_NI = PSI_NI + QOPI*(PSNI*SSUM + PDNI*SDIF);
% 
%          QTANM = QTANM + QOPI*(PSNI*SSUM + PDNI*SDIF);
% 
%          DQDM(JM) = DQDM(JM) + QOPI*(-PSNI*DSIM + PDNI*DSIM);
%          DQDM(JO) = DQDM(JO) + QOPI*(-PSNI*DSIO - PDNI*DSIO);
%          DQDM(JP) = DQDM(JP) + QOPI*( PSNI*(DSIO+DSIM)            + PDNI*(DSIO-DSIM));
% % C
% % C
% % C------- calculate source contribution to Psi  for  0-2  half-panel
%          DXINV = 1.0/(X0-X2);
%          PSUM = X2*(T2-APAN) - X0*(T0-APAN) + 0.5*YY*(G0-G2);
%          PDIF = ((X0+X2)*PSUM + RS0*(T0-APAN) - RS2*(T2-APAN)     + (X2-X0)*YY) * DXINV;
% 
%          PSX0 =  -(T0-APAN);
%          PSX2 =    T2-APAN;
%          PSYY =  0.5*(G0-G2);
% 
%          PDX0 = ((X0+X2)*PSX0 + PSUM + 2.0*X0*(T0-APAN) - PDIF) * DXINV;
%          PDX2 = ((X0+X2)*PSX2 + PSUM - 2.0*X2*(T2-APAN) + PDIF) * DXINV;
%          PDYY = ((X0+X2)*PSYY + 2.0*(X2-X0 + YY*(T0-T2))      ) * DXINV;
% 
%          DSP = sqrt((X(JQ)-X(JO))^2 + (Y(JQ)-Y(JO))^2);
%          DSIP = 1.0/DSP;
% 
% % CCC         SIG2 = (SIG(JQ) - SIG(JO))*DSIP
% % CCC         SIG0 = (SIG(JP) - SIG(JO))*DSIO
% % CCC         SSUM = SIG2 + SIG0
% % CCC         SDIF = SIG2 - SIG0
% % C
%          SSUM = (SIG(JQ) - SIG(JO))*DSIP + (SIG(JP) - SIG(JO))*DSIO;
%          SDIF = (SIG(JQ) - SIG(JO))*DSIP - (SIG(JP) - SIG(JO))*DSIO;
% 
%          PSI = PSI + QOPI*(PSUM*SSUM + PDIF*SDIF);
% 
% %C------- dPsi/dm
%          DZDM(JO) = DZDM(JO) + QOPI*(-PSUM*(DSIP+DSIO)                          - PDIF*(DSIP-DSIO));
%          DZDM(JP) = DZDM(JP) + QOPI*( PSUM*DSIO - PDIF*DSIO);
%          DZDM(JQ) = DZDM(JQ) + QOPI*( PSUM*DSIP + PDIF*DSIP);
% 
% %C------- dPsi/dni
%          PSNI = PSX0*(X1I+X2I)*0.5 + PSX2*X2I + PSYY*YYI;
%          PDNI = PDX0*(X1I+X2I)*0.5 + PDX2*X2I + PDYY*YYI;
%          PSI_NI = PSI_NI + QOPI*(PSNI*SSUM + PDNI*SDIF);
% 
%          QTANM = QTANM + QOPI*(PSNI*SSUM + PDNI*SDIF);
% 
%          DQDM(JO) = DQDM(JO) + QOPI*(-PSNI*(DSIP+DSIO)                      - PDNI*(DSIP-DSIO));
%          DQDM(JP) = DQDM(JP) + QOPI*( PSNI*DSIO - PDNI*DSIO);
%          DQDM(JQ) = DQDM(JQ) + QOPI*( PSNI*DSIP + PDNI*DSIP);
% 
%         end
% % C
% % C------ calculate vortex panel contribution to Psi
%         DXINV = 1.0/(X1-X2);
%         PSIS = 0.5*X1*G1 - 0.5*X2*G2 + X2 - X1 + YY*(T1-T2);
%         PSID = ((X1+X2)*PSIS + 0.5*(RS2*G2-RS1*G1 + X1*X1-X2*X2))*DXINV;
% 
%         PSX1 = 0.5*G1;
%         PSX2 = -.5*G2;
%         PSYY = T1-T2;
% 
%         PDX1 = ((X1+X2)*PSX1 + PSIS - X1*G1 - PSID)*DXINV;
%         PDX2 = ((X1+X2)*PSX2 + PSIS + X2*G2 + PSID)*DXINV;
%         PDYY = ((X1+X2)*PSYY - YY*(G1-G2)         )*DXINV;
% 
%         GSUM1 = GAMU(JP,1) + GAMU(JO,1);
%         GSUM2 = GAMU(JP,2) + GAMU(JO,2);
%         GDIF1 = GAMU(JP,1) - GAMU(JO,1);
%         GDIF2 = GAMU(JP,2) - GAMU(JO,2);
% 
%         GSUM = GAM(JP) + GAM(JO);
%         GDIF = GAM(JP) - GAM(JO);
% 
%         PSI = PSI + QOPI*(PSIS*GSUM + PSID*GDIF);
% 
% %C------ dPsi/dGam
%         DZDG(JO) = DZDG(JO) + QOPI*(PSIS-PSID);
%         DZDG(JP) = DZDG(JP) + QOPI*(PSIS+PSID);
% 
% %C------ dPsi/dni
%         PSNI = PSX1*X1I + PSX2*X2I + PSYY*YYI;
%         PDNI = PDX1*X1I + PDX2*X2I + PDYY*YYI;
%         PSI_NI = PSI_NI + QOPI*(GSUM*PSNI + GDIF*PDNI);
% 
%         QTAN1 = QTAN1 + QOPI*(GSUM1*PSNI + GDIF1*PDNI);
%         QTAN2 = QTAN2 + QOPI*(GSUM2*PSNI + GDIF2*PDNI);
% 
%         DQDG(JO) = DQDG(JO) + QOPI*(PSNI - PDNI);
%         DQDG(JP) = DQDG(JP) + QOPI*(PSNI + PDNI);
% 
%         if(GEOLIN==true) 
% % C
% % C------- dPsi/dn
%          DZDN(JO) = DZDN(JO)+ QOPI*GSUM*(PSX1*X1O + PSX2*X2O + PSYY*YYO)      + QOPI*GDIF*(PDX1*X1O + PDX2*X2O + PDYY*YYO);
%          DZDN(JP) = DZDN(JP)+ QOPI*GSUM*(PSX1*X1P + PSX2*X2P + PSYY*YYP)     + QOPI*GDIF*(PDX1*X1P + PDX2*X2P + PDYY*YYP);
% %C------- dPsi/dP
%          Z_QDOF0 = Z_QDOF0      + QOPI*((PSIS-PSID)*QF0(JO) + (PSIS+PSID)*QF0(JP));
%          Z_QDOF1 = Z_QDOF1      + QOPI*((PSIS-PSID)*QF1(JO) + (PSIS+PSID)*QF1(JP));
%          Z_QDOF2 = Z_QDOF2     + QOPI*((PSIS-PSID)*QF2(JO) + (PSIS+PSID)*QF2(JP));
%          Z_QDOF3 = Z_QDOF3        + QOPI*((PSIS-PSID)*QF3(JO) + (PSIS+PSID)*QF3(JP));
%         end
%         
%       end
% 
%       PSIG = 0.5*YY*(G1-G2) + X2*(T2-APAN) - X1*(T1-APAN);
%       PGAM = 0.5*X1*G1 - 0.5*X2*G2 + X2 - X1 + YY*(T1-T2);
% 
%       PSIGX1 = -(T1-APAN);
%       PSIGX2 =   T2-APAN;
%       PSIGYY = 0.5*(G1-G2);
%       PGAMX1 = 0.5*G1;
%       PGAMX2 = -.5*G2;
%       PGAMYY = T1-T2;
% 
%       PSIGNI = PSIGX1*X1I + PSIGX2*X2I + PSIGYY*YYI;
%       PGAMNI = PGAMX1*X1I + PGAMX2*X2I + PGAMYY*YYI;
% % C
% % C---- TE panel source and vortex strengths
%       SIGTE1 = 0.5*SCS*(GAMU(JP,1) - GAMU(JO,1));
%       SIGTE2 = 0.5*SCS*(GAMU(JP,2) - GAMU(JO,2));
%       GAMTE1 = -.5*SDS*(GAMU(JP,1) - GAMU(JO,1));
%       GAMTE2 = -.5*SDS*(GAMU(JP,2) - GAMU(JO,2));
% 
%       SIGTE = 0.5*SCS*(GAM(JP) - GAM(JO));
%       GAMTE = -.5*SDS*(GAM(JP) - GAM(JO));
% 
% %C---- TE panel contribution to Psi
%       PSI = PSI + HOPI*(PSIG*SIGTE + PGAM*GAMTE);
% 
% %C---- dPsi/dGam
%       DZDG(JO) = DZDG(JO) - HOPI*PSIG*SCS*0.5;
%       DZDG(JP) = DZDG(JP) + HOPI*PSIG*SCS*0.5;
% 
%       DZDG(JO) = DZDG(JO) + HOPI*PGAM*SDS*0.5;
%       DZDG(JP) = DZDG(JP) - HOPI*PGAM*SDS*0.5;
% 
% %C---- dPsi/dni
%       PSI_NI = PSI_NI + HOPI*(PSIGNI*SIGTE + PGAMNI*GAMTE);
% 
%       QTAN1 = QTAN1 + HOPI*(PSIGNI*SIGTE1 + PGAMNI*GAMTE1);
%       QTAN2 = QTAN2 + HOPI*(PSIGNI*SIGTE2 + PGAMNI*GAMTE2);
% 
%       DQDG(JO) = DQDG(JO) - HOPI*(PSIGNI*0.5*SCS - PGAMNI*0.5*SDS);
%       DQDG(JP) = DQDG(JP) + HOPI*(PSIGNI*0.5*SCS - PGAMNI*0.5*SDS);
% 
%       if(GEOLIN==true) 
% 
% %C----- dPsi/dn
%        DZDN(JO) = DZDN(JO)+ HOPI*(PSIGX1*X1O + PSIGX2*X2O + PSIGYY*YYO)*SIGTE    + HOPI*(PGAMX1*X1O + PGAMX2*X2O + PGAMYY*YYO)*GAMTE;
%        DZDN(JP) = DZDN(JP)       + HOPI*(PSIGX1*X1P + PSIGX2*X2P + PSIGYY*YYP)*SIGTE      + HOPI*(PGAMX1*X1P + PGAMX2*X2P + PGAMYY*YYP)*GAMTE;
% % C
% % C----- dPsi/dP
%        Z_QDOF0 = Z_QDOF0 + HOPI*PSIG*0.5*(QF0(JP)-QF0(JO))*SCS                 - HOPI*PGAM*0.5*(QF0(JP)-QF0(JO))*SDS;
%        Z_QDOF1 = Z_QDOF1 + HOPI*PSIG*0.5*(QF1(JP)-QF1(JO))*SCS              - HOPI*PGAM*0.5*(QF1(JP)-QF1(JO))*SDS;
%        Z_QDOF2 = Z_QDOF2 + HOPI*PSIG*0.5*(QF2(JP)-QF2(JO))*SCS                   - HOPI*PGAM*0.5*(QF2(JP)-QF2(JO))*SDS;
%        Z_QDOF3 = Z_QDOF3 + HOPI*PSIG*0.5*(QF3(JP)-QF3(JO))*SCS                 - HOPI*PGAM*0.5*(QF3(JP)-QF3(JO))*SDS;
%       end
% 
%              
%                    PSI = PSI + QINF*(COSA*YI - SINA*XI);
% % C
% % C---- dPsi/dn
%       PSI_NI = PSI_NI + QINF*(COSA*NYI - SINA*NXI);
% 
%       QTAN1 = QTAN1 + QINF*NYI;
%       QTAN2 = QTAN2 - QINF*NXI;
% % C
% % C---- dPsi/dQinf
%       Z_QINF = Z_QINF + (COSA*YI - SINA*XI);
% % C
% % C---- dPsi/dalfa
%       Z_ALFA = Z_ALFA - QINF*(SINA*YI + COSA*XI);
% % C
%       if (LIMAGE==false) 
%       return
%       end
%       
%       for JO=1: N
%         JP = JO+1;
% 
%         JM = JO-1;
%         JQ = JP+1;
% 
%         if(JO==1)
%          JM = JO;
%         elseif (JO==N-1) 
%          JQ = JP;
%         elseif (JO==N) 
%          JP = 1;
%          if((X(JO)-X(JP))^2 + (Y(JO)-Y(JP))^2 < SEPS^2) 
%            return
%          end
%         end
% 
%         DSO = sqrt((X(JO)-X(JP))^2 + (Y(JO)-Y(JP))^2)
% % C
% % C------ skip null panel
%         if (DSO == 0.0) 
%             continue
%         end
% 
%         DSIO = 1.0 / DSO;
% 
%         APAN = pi - APANEL(JO) + 2.0*ALFA;
% 
%         XJO = X(JO) + 2.0*(YIMAGE+Y(JO))*SINA;
%         YJO = Y(JO) - 2.0*(YIMAGE+Y(JO))*COSA;
%         XJP = X(JP) + 2.0*(YIMAGE+Y(JP))*SINA;
%         YJP = Y(JP) - 2.0*(YIMAGE+Y(JP))*COSA;
% 
%         RX1 = XI - XJO;
%         RY1 = YI - YJO;
%         RX2 = XI - XJP;
%         RY2 = YI - YJP;
% 
%         SX = (XJP - XJO) * DSIO;
%         SY = (YJP - YJO) * DSIO;
% 
%         X1 = SX*RX1 + SY*RY1;
%         X2 = SX*RX2 + SY*RY2;
%         YY = SX*RY1 - SY*RX1;
% 
%         RS1 = RX1*RX1 + RY1*RY1;
%         RS2 = RX2*RX2 + RY2*RY2;
% % C
% % C------ set reflection flag SGN to avoid branch problems with arctan
%         if(IO>=1 && IO~=N) 
% %C------- no problem on airfoil surface
%          SGN = 1.0;
%         else
% % C------- make sure arctan falls between  -/+  Pi/2
%          SGN = 1*sign(YY);
%         end
% 
% % C------ set log(r^2) and arctan(p.x/p.y), correcting for reflection if any
%         G1 = log(RS1);
%         T1 = atan2(SGN*X1,SGN*YY) + (0.5 - 0.5*SGN)*pi;
% 
%         G2 = log(RS2);
%         T2 = atan2(SGN*X2,SGN*YY) + (0.5 - 0.5*SGN)*pi
% 
%         X1I = SX*NXI + SY*NYI;
%         X2I = SX*NXI + SY*NYI;
%         YYI = SX*NYI - SY*NXI;
% 
%         if (GEOLIN==true)
%          NXO = NX(JO);
%          NYO = NY(JO);
%          NXP = NX(JP);
%          NYP = NY(JP);
% 
%          X1O =-((RX1-X1*SX)*NXO + (RY1-X1*SY)*NYO)*DSIO-(SX*NXO+SY*NYO);
%          X1P = ((RX1-X1*SX)*NXP + (RY1-X1*SY)*NYP)*DSIO;
%          X2O =-((RX2-X2*SX)*NXO + (RY2-X2*SY)*NYO)*DSIO;
%          X2P = ((RX2-X2*SX)*NXP + (RY2-X2*SY)*NYP)*DSIO-(SX*NXP+SY*NYP);
%          YYO = ((RX1+X1*SY)*NYO - (RY1-X1*SX)*NXO)*DSIO-(SX*NYO-SY*NXO);
%          YYP =-((RX1-X1*SY)*NYP - (RY1+X1*SX)*NXP)*DSIO;
%         end
% 
%         if (JO==N) 
%             break
%         end
% 
%         if (SIGLIN==true) 
% % C
% % C------- set up midpoint quantities
%          X0 = 0.5*(X1+X2);
%          RS0 = X0*X0 + YY*YY;
%          G0 = log(RS0);
%          T0 = atan2(SGN*X0,SGN*YY) + (0.5 - 0.5*SGN)*pi;
% % C
% % C------- calculate source contribution to Psi  for  1-0  half-panel
%          DXINV = 1.0/(X1-X0);
%          PSUM = X0*(T0-APAN) - X1*(T1-APAN) + 0.5*YY*(G1-G0);
%          PDIF = ((X1+X0)*PSUM + RS1*(T1-APAN) - RS0*(T0-APAN)    + (X0-X1)*YY) * DXINV;
% 
%          PSX1 =  -(T1-APAN);
%          PSX0 =    T0-APAN;
%          PSYY =  0.5*(G1-G0);
% 
%          PDX1 = ((X1+X0)*PSX1 + PSUM + 2.0*X1*(T1-APAN) - PDIF) * DXINV;
%          PDX0 = ((X1+X0)*PSX0 + PSUM - 2.0*X0*(T0-APAN) + PDIF) * DXINV;
%          PDYY = ((X1+X0)*PSYY + 2.0*(X0-X1 + YY*(T1-T0))      ) * DXINV;
% 
%          DSM = sqrt((X(JP)-X(JM))^2 + (Y(JP)-Y(JM))^2);
%          DSIM = 1.0/DSM;
% 
% % CCC      SIG0 = (SIG(JP) - SIG(JO))*DSIO
% % CCC      SIG1 = (SIG(JP) - SIG(JM))*DSIM
% % CCC      SSUM = SIG0 + SIG1
% % CCC      SDIF = SIG0 - SIG1
% % C
%          SSUM = (SIG(JP) - SIG(JO))*DSIO + (SIG(JP) - SIG(JM))*DSIM;
%          SDIF = (SIG(JP) - SIG(JO))*DSIO - (SIG(JP) - SIG(JM))*DSIM;
% 
%          PSI = PSI + QOPI*(PSUM*SSUM + PDIF*SDIF);
% % C
% % C------- dPsi/dm
%          DZDM(JM) = DZDM(JM) + QOPI*(-PSUM*DSIM + PDIF*DSIM);
%          DZDM(JO) = DZDM(JO) + QOPI*(-PSUM*DSIO - PDIF*DSIO);
%          DZDM(JP) = DZDM(JP) + QOPI*( PSUM*(DSIO+DSIM)     + PDIF*(DSIO-DSIM));
% % C
% % C------- dPsi/dni
%          PSNI = PSX1*X1I + PSX0*(X1I+X2I)*0.5 + PSYY*YYI;
%          PDNI = PDX1*X1I + PDX0*(X1I+X2I)*0.5 + PDYY*YYI;
%          PSI_NI = PSI_NI + QOPI*(PSNI*SSUM + PDNI*SDIF);
% 
%          QTANM = QTANM + QOPI*(PSNI*SSUM + PDNI*SDIF);
% 
%          DQDM(JM) = DQDM(JM) + QOPI*(-PSNI*DSIM + PDNI*DSIM);
%          DQDM(JO) = DQDM(JO) + QOPI*(-PSNI*DSIO - PDNI*DSIO);
%          DQDM(JP) = DQDM(JP) + QOPI*( PSNI*(DSIO+DSIM) + PDNI*(DSIO-DSIM));
% 
% 
% %C------- calculate source contribution to Psi  for  0-2  half-panel
%          DXINV = 1.0/(X0-X2);
%          PSUM = X2*(T2-APAN) - X0*(T0-APAN) + 0.5*YY*(G0-G2);
%          PDIF = ((X0+X2)*PSUM + RS0*(T0-APAN) - RS2*(T2-APAN) + (X2-X0)*YY) * DXINV;
% 
%          PSX0 =  -(T0-APAN);
%          PSX2 =    T2-APAN;
%          PSYY =  0.5*(G0-G2);
% 
%          PDX0 = ((X0+X2)*PSX0 + PSUM + 2.0*X0*(T0-APAN) - PDIF) * DXINV;
%          PDX2 = ((X0+X2)*PSX2 + PSUM - 2.0*X2*(T2-APAN) + PDIF) * DXINV;
%          PDYY = ((X0+X2)*PSYY + 2.0*(X2-X0 + YY*(T0-T2))      ) * DXINV;
% 
%          DSP = sqrt((X(JQ)-X(JO))^2 + (Y(JQ)-Y(JO))^2);
%          DSIP = 1.0/DSP;
% % 
% % CCC         SIG2 = (SIG(JQ) - SIG(JO))*DSIP
% % CCC         SIG0 = (SIG(JP) - SIG(JO))*DSIO
% % CCC         SSUM = SIG2 + SIG0
% % CCC         SDIF = SIG2 - SIG0
% % C
%          SSUM = (SIG(JQ) - SIG(JO))*DSIP + (SIG(JP) - SIG(JO))*DSIO;
%          SDIF = (SIG(JQ) - SIG(JO))*DSIP - (SIG(JP) - SIG(JO))*DSIO;
% 
%          PSI = PSI + QOPI*(PSUM*SSUM + PDIF*SDIF);
% 
% %C------- dPsi/dm
%          DZDM(JO) = DZDM(JO) + QOPI*(-PSUM*(DSIP+DSIO)      - PDIF*(DSIP-DSIO));
%          DZDM(JP) = DZDM(JP) + QOPI*( PSUM*DSIO - PDIF*DSIO);
%          DZDM(JQ) = DZDM(JQ) + QOPI*( PSUM*DSIP + PDIF*DSIP);
% 
% %C------- dPsi/dni
%          PSNI = PSX0*(X1I+X2I)*0.5 + PSX2*X2I + PSYY*YYI;
%          PDNI = PDX0*(X1I+X2I)*0.5 + PDX2*X2I + PDYY*YYI;
%          PSI_NI = PSI_NI + QOPI*(PSNI*SSUM + PDNI*SDIF);
% 
%          QTANM = QTANM + QOPI*(PSNI*SSUM + PDNI*SDIF);
% 
%          DQDM(JO) = DQDM(JO) + QOPI*(-PSNI*(DSIP+DSIO)          - PDNI*(DSIP-DSIO));
%          DQDM(JP) = DQDM(JP) + QOPI*( PSNI*DSIO - PDNI*DSIO);
%          DQDM(JQ) = DQDM(JQ) + QOPI*( PSNI*DSIP + PDNI*DSIP);
% 
%         end
% 
% %C------ calculate vortex panel contribution to Psi
%         DXINV = 1.0/(X1-X2);
%         PSIS = 0.5*X1*G1 - 0.5*X2*G2 + X2 - X1 + YY*(T1-T2);
%         PSID = ((X1+X2)*PSIS + 0.5*(RS2*G2-RS1*G1 + X1*X1-X2*X2))*DXINV;
% 
%         PSX1 = 0.5*G1;
%         PSX2 = -.5*G2;
%         PSYY = T1-T2;
% 
%         PDX1 = ((X1+X2)*PSX1 + PSIS - X1*G1 - PSID)*DXINV;
%         PDX2 = ((X1+X2)*PSX2 + PSIS + X2*G2 + PSID)*DXINV;
%         PDYY = ((X1+X2)*PSYY - YY*(G1-G2)         )*DXINV;
% 
%         GSUM1 = GAMU(JP,1) + GAMU(JO,1);
%         GSUM2 = GAMU(JP,2) + GAMU(JO,2);
%         GDIF1 = GAMU(JP,1) - GAMU(JO,1);
%         GDIF2 = GAMU(JP,2) - GAMU(JO,2);
% 
%         GSUM = GAM(JP) + GAM(JO);
%         GDIF = GAM(JP) - GAM(JO);
% 
%         PSI = PSI - QOPI*(PSIS*GSUM + PSID*GDIF);
% 
% %C------ dPsi/dGam
%         DZDG(JO) = DZDG(JO) - QOPI*(PSIS-PSID);
%         DZDG(JP) = DZDG(JP) - QOPI*(PSIS+PSID);
% 
% %C------ dPsi/dni
%         PSNI = PSX1*X1I + PSX2*X2I + PSYY*YYI;
%         PDNI = PDX1*X1I + PDX2*X2I + PDYY*YYI;
%         PSI_NI = PSI_NI - QOPI*(GSUM*PSNI + GDIF*PDNI);
% 
%         QTAN1 = QTAN1 - QOPI*(GSUM1*PSNI + GDIF1*PDNI);
%         QTAN2 = QTAN2 - QOPI*(GSUM2*PSNI + GDIF2*PDNI);
% 
%         DQDG(JO) = DQDG(JO) - QOPI*(PSNI - PDNI);
%         DQDG(JP) = DQDG(JP) - QOPI*(PSNI + PDNI);
% 
%         if (GEOLIN==true) 
% % C
% % C------- dPsi/dn
%          DZDN(JO) = DZDN(JO)- QOPI*GSUM*(PSX1*X1O + PSX2*X2O + PSYY*YYO) - QOPI*GDIF*(PDX1*X1O + PDX2*X2O + PDYY*YYO);
%          DZDN(JP) = DZDN(JP)- QOPI*GSUM*(PSX1*X1P + PSX2*X2P + PSYY*YYP) - QOPI*GDIF*(PDX1*X1P + PDX2*X2P + PDYY*YYP);
% %C------- dPsi/dP
%          Z_QDOF0 = Z_QDOF0         - QOPI*((PSIS-PSID)*QF0(JO) + (PSIS+PSID)*QF0(JP));
%          Z_QDOF1 = Z_QDOF1      - QOPI*((PSIS-PSID)*QF1(JO) + (PSIS+PSID)*QF1(JP));
%          Z_QDOF2 = Z_QDOF2       - QOPI*((PSIS-PSID)*QF2(JO) + (PSIS+PSID)*QF2(JP));
%          Z_QDOF3 = Z_QDOF3       - QOPI*((PSIS-PSID)*QF3(JO) + (PSIS+PSID)*QF3(JP));
%         end
% 
% 
%       end
%          
%       PSIG = 0.5*YY*(G1-G2) + X2*(T2-APAN) - X1*(T1-APAN);
%       PGAM = 0.5*X1*G1 - 0.5*X2*G2 + X2 - X1 + YY*(T1-T2);
% 
%       PSIGX1 = -(T1-APAN);
%       PSIGX2 =   T2-APAN;
%       PSIGYY = 0.5*(G1-G2);
%       PGAMX1 = 0.5*G1;
%       PGAMX2 = -.5*G2;
%       PGAMYY = T1-T2;
% 
%       PSIGNI = PSIGX1*X1I + PSIGX2*X2I + PSIGYY*YYI;
%       PGAMNI = PGAMX1*X1I + PGAMX2*X2I + PGAMYY*YYI;
% 
% %C---- TE panel source and vortex strengths
%       SIGTE1 = 0.5*SCS*(GAMU(JP,1) - GAMU(JO,1));
%       SIGTE2 = 0.5*SCS*(GAMU(JP,2) - GAMU(JO,2));
%       GAMTE1 = -.5*SDS*(GAMU(JP,1) - GAMU(JO,1));
%       GAMTE2 = -.5*SDS*(GAMU(JP,2) - GAMU(JO,2));
% 
%       SIGTE = 0.5*SCS*(GAM(JP) - GAM(JO));
%       GAMTE = -.5*SDS*(GAM(JP) - GAM(JO));
% 
% %C---- TE panel contribution to Psi
%       PSI = PSI + HOPI*(PSIG*SIGTE - PGAM*GAMTE);
% 
% %C---- dPsi/dGam
%       DZDG(JO) = DZDG(JO) - HOPI*PSIG*SCS*0.5;
%       DZDG(JP) = DZDG(JP) + HOPI*PSIG*SCS*0.5;
% 
%       DZDG(JO) = DZDG(JO) - HOPI*PGAM*SDS*0.5;
%       DZDG(JP) = DZDG(JP) + HOPI*PGAM*SDS*0.5;
% 
% %C---- dPsi/dni
%       PSI_NI = PSI_NI + HOPI*(PSIGNI*SIGTE - PGAMNI*GAMTE);
% 
%       QTAN1 = QTAN1 + HOPI*(PSIGNI*SIGTE1 - PGAMNI*GAMTE1);
%       QTAN2 = QTAN2 + HOPI*(PSIGNI*SIGTE2 - PGAMNI*GAMTE2);
% 
%       DQDG(JO) = DQDG(JO) - HOPI*(PSIGNI*0.5*SCS + PGAMNI*0.5*SDS);
%       DQDG(JP) = DQDG(JP) + HOPI*(PSIGNI*0.5*SCS + PGAMNI*0.5*SDS);
% 
%       if (GEOLIN==true) 
% 
% %C----- dPsi/dn
%        DZDN(JO) = DZDN(JO)  + HOPI*(PSIGX1*X1O + PSIGX2*X2O + PSIGYY*YYO)*SIGTE  - HOPI*(PGAMX1*X1O + PGAMX2*X2O + PGAMYY*YYO)*GAMTE;
%        DZDN(JP) = DZDN(JP) + HOPI*(PSIGX1*X1P + PSIGX2*X2P + PSIGYY*YYP)*SIGTE      - HOPI*(PGAMX1*X1P + PGAMX2*X2P + PGAMYY*YYP)*GAMTE;
% % C
% % C----- dPsi/dP
%        Z_QDOF0 = Z_QDOF0 + HOPI*PSIG*0.5*(QF0(JP)-QF0(JO))*SCS                + HOPI*PGAM*0.5*(QF0(JP)-QF0(JO))*SDS;
%        Z_QDOF1 = Z_QDOF1 + HOPI*PSIG*0.5*(QF1(JP)-QF1(JO))*SCS                + HOPI*PGAM*0.5*(QF1(JP)-QF1(JO))*SDS;
%        Z_QDOF2 = Z_QDOF2 + HOPI*PSIG*0.5*(QF2(JP)-QF2(JO))*SCS                + HOPI*PGAM*0.5*(QF2(JP)-QF2(JO))*SDS;
%        Z_QDOF3 = Z_QDOF3 + HOPI*PSIG*0.5*(QF3(JP)-QF3(JO))*SCS               + HOPI*PGAM*0.5*(QF3(JP)-QF3(JO))*SDS;
%       end
%          
% return
% end
% C-----------------------------------------------------------------------
% C     Calculates current streamfunction Psi at panel node or wake node
% C     I due to freestream and all bound vorticity Gam on the airfoil. 
% C     Sensitivities of Psi with respect to alpha (Z_ALFA) and inverse
% C     Qspec DOFs (Z_QDOF0,Z_QDOF1) which influence Gam in inverse cases.
% C     Also calculates the sensitivity vector dPsi/dGam (DZDG).
% C
% C     If SIGLIN=True, then Psi includes the effects of the viscous
% C     source distribution Sig and the sensitivity vector dPsi/dSig
% C     (DZDM) is calculated.
% C
% C     If GEOLIN=True, then the geometric sensitivity vector dPsi/dn
% C     is calculated, where p.n is the normal motion of the jth node.
% C
% C          Airfoil:  1   < I < N
% C          Wake:     N+1 < I < N+NW
% C-----------------------------------------------------------------------
% C
% C---- distance tolerance for determining if two pos are the same

 seps = (p.s(p.n)-p.s(1)) * 0.00001;
	qopi = 0.25/pi;
    hopi = 0.50/pi;
 qinf = 1.0;
 limage=false;

	io = i;
    

     
    
    
	
	cosa = cos(p.alfa);
	sina = sin(p.alfa);


	 jp = 0;%// added arcds
	 apan = 0.0;
	 yy = 0.0;
	 g1 = 0.0;
	 g2 = 0.0;
	 x1i = 0.0;
	 x2i = 0.0;
	 yyi = 0.0;
	 x1o = 0.0;
	 x1p = 0.0;
	 x2o = 0.0;
	 x2p = 0.0;
	 yyo = 0.0;
	 yyp = 0.0;




	for (jo=1:p.n)
		p.dzdg(jo) = 0.0;
		p.dzdn(jo) = 0.0;
		p.dqdg(jo) = 0.0;
    end
	
	for (jo=1:p.n)
		p.dzdm(jo) = 0.0;
		p.dqdm(jo) = 0.0;
    end
	
	p.z_qinf = 0.0;
	p.z_alfa = 0.0;
	p.z_qdof0 = 0.0;
	p.z_qdof1 = 0.0;
	p.z_qdof2 = 0.0;
	p.z_qdof3 = 0.0;
	
	psi  = 0.0;
	psi_ni = 0.0;
	
	p.qtan1 = 0.0;
	p.qtan2 = 0.0;
	p.qtanm = 0.0;

	
	if(p.sharp)
		scs = 1.0;
		sds = 0.0;
    
    else
		scs = p.ante/p.dste;
		sds = p.aste/p.dste;
    end
	for (jo=1:p.n)
		jp = jo+1;
		jm = jo-1;
		jq = jp+1;
		
		if(jo==1) 
            jm = jo;
		else 
			if(jo==p.n-1) jq = jp;
            else
				if(jo==p.n) 
					jp = 1;
					if((p.x(jo)-p.x(jp))*(p.x(jo)-p.x(jp)) + (p.y(jo)-p.y(jp))*(p.y(jo)-p.y(jp))< seps*seps) 
						%goto stop12;
                        
                        
                        
                        
                        	%//**** freestream terms
	psi = psi + qinf*(cosa*yi - sina*xi);
	
	%//---- dpsi/dn
	psi_ni = psi_ni + qinf*(cosa*nyi - sina*nxi);
	
	p.qtan1 = p.qtan1 + qinf*nyi;
	p.qtan2 = p.qtan2 - qinf*nxi;
	
	%//---- dpsi/dqinf
	p.z_qinf = p.z_qinf + (cosa*yi - sina*xi);
	
	%//---- dpsi/dalfa
	p.z_alfa = p.z_alfa - qinf*(sina*yi + cosa*xi);
	
	if(limage==false)
		return;
    end

	for(jo=1:p.n)
		jp = jo+1;
		
		jm = jo-1;
		jq = jp+1;
		
		if(jo==1) 
            jm = jo;
        else
			if(jo==p.n-1)  jq = jp;
			else 
				if(jo==p.n) jp = 1;
				if((p.x(jo)-p.x(jp))*(p.x(jo)-p.x(jp)) + (p.y(jo)-p.y(jp))*(p.y(jo)-p.y(jp)) < seps*seps)
					return;
                end
                end
            end
        end
		
		dso = sqrt((p.x(jo)-p.x(jp))*(p.x(jo)-p.x(jp)) + (p.y(jo)-p.y(jp))*(p.y(jo)-p.y(jp)));
		
		%//------ skip null panel
		if (dso == 0.0) 
            %goto stop20; %// check - unsafe condition
            
            		nothing = 1;
                    continue;
        end

            
    
            
            
            
            
		dsio = 1.0 / dso;
		%//%//ccc   apan = p.apanel(jo)
		apan = pi - p.apanel(jo) + 2.0*p.alfa;
		
		xjo = p.x(jo) + 2.0*(yimage+p.y(jo))*sina;
		yjo = p.y(jo) - 2.0*(yimage+p.y(jo))*cosa;
		xjp = p.x(jp) + 2.0*(yimage+p.y(jp))*sina;
		yjp = p.y(jp) - 2.0*(yimage+p.y(jp))*cosa;
		
		rx1 = xi - xjo;
		ry1 = yi - yjo;
		rx2 = xi - xjp;
		ry2 = yi - yjp;
		
		sx = (xjp - xjo) * dsio;
		sy = (yjp - yjo) * dsio;
		
		x1 = sx*rx1 + sy*ry1;
		x2 = sx*rx2 + sy*ry2;
		yy = sx*ry1 - sy*rx1;
		
		rs1 = rx1*rx1 + ry1*ry1;
		rs2 = rx2*rx2 + ry2*ry2;
		
		%//------ set reflection flag sgn to avoid branch problems with arctan
		if(io>=1 && io<=p.n) 
			%//------- no problem on airfoil surface
			sgn = 1.0;
		
        else
			%//------- make sure arctan falls between	-/+  pi/2
			sgn = 1*sign(yy);
        end
		
		%//------ set log(r^2) and arctan(p.x/p.y), correcting for reflection if any
		g1 = log(rs1);
		t1 = atan2(sgn*x1,sgn*yy) + (0.5 - 0.5*sgn)*pi;
		
		g2 = log(rs2);
		t2 = atan2(sgn*x2,sgn*yy) + (0.5 - 0.5*sgn)*pi;
		
		x1i = sx*nxi + sy*nyi;
		x2i = sx*nxi + sy*nyi;
		yyi = sx*nyi - sy*nxi;
		
		if(geolin) 
			nxo = p.nx(jo);
			nyo = p.ny(jo);
			nxp = p.nx(jp);
			nyp = p.ny(jp);
			
			x1o =-((rx1-x1*sx)*nxo + (ry1-x1*sy)*nyo)*dsio-(sx*nxo+sy*nyo);
			x1p = ((rx1-x1*sx)*nxp + (ry1-x1*sy)*nyp)*dsio;
			x2o =-((rx2-x2*sx)*nxo + (ry2-x2*sy)*nyo)*dsio;
			x2p = ((rx2-x2*sx)*nxp + (ry2-x2*sy)*nyp)*dsio-(sx*nxp+sy*nyp);
			yyo = ((rx1+x1*sy)*nyo - (ry1-x1*sx)*nxo)*dsio-(sx*nyo-sy*nxo);
			yyp =-((rx1-x1*sy)*nyp - (ry1+x1*sx)*nxp)*dsio;
        end
		
		if(jo==p.n) 
            break;
        end
			
		if(siglin) 
			
			%//------- set up midpo quantities
			x0 = 0.5*(x1+x2);
			rs0 = x0*x0 + yy*yy;
			g0 = log(rs0);
			t0 = atan2(sgn*x0,sgn*yy) + (0.5 - 0.5*sgn)*pi;
			
			%//------- calculate source contribution to psi	for  1-0  half-panel
			dxinv = 1.0/(x1-x0);
			psum = x0*(t0-apan) - x1*(t1-apan) + 0.5*yy*(g1-g0);
			pdif = ((x1+x0)*psum + rs1*(t1-apan) - rs0*(t0-apan)				+ (x0-x1)*yy) * dxinv;
			
			psx1 =  -(t1-apan);
			psx0 =	t0-apan;
			psyy =  0.5*(g1-g0);
			
			pdx1 = ((x1+x0)*psx1 + psum + 2.0*x1*(t1-apan) - pdif) * dxinv;
			pdx0 = ((x1+x0)*psx0 + psum - 2.0*x0*(t0-apan) + pdif) * dxinv;
			pdyy = ((x1+x0)*psyy + 2.0*(x0-x1 + yy*(t1-t0))	  ) * dxinv;
			
			dsm = sqrt((p.x(jp)-p.x(jm))*(p.x(jp)-p.x(jm)) + (p.y(jp)-p.y(jm))*(p.y(jp)-p.y(jm)));
			dsim = 1.0/dsm;
			
			%//%//ccc    sig0 = (p.sig(jp) - p.sig(jo))*dsio
			%//%//ccc    sig1 = (p.sig(jp) - p.sig(jm))*dsim
			%//%//ccc    ssum = sig0 + sig1
			%//%//ccc    sdif = sig0 - sig1
			
			ssum = (p.sig(jp) - p.sig(jo))*dsio + (p.sig(jp) - p.sig(jm))*dsim;
			sdif = (p.sig(jp) - p.sig(jo))*dsio - (p.sig(jp) - p.sig(jm))*dsim;
			
			psi = psi + qopi*(psum*ssum + pdif*sdif);
			
			%//------- dpsi/dm
			p.dzdm(jm) = p.dzdm(jm) + qopi*(-psum*dsim + pdif*dsim);
			p.dzdm(jo) = p.dzdm(jo) + qopi*(-psum*dsio - pdif*dsio);
			p.dzdm(jp) = p.dzdm(jp) + qopi*( psum*(dsio+dsim)+ pdif*(dsio-dsim));
			
			%//------- dpsi/dni
			psni = psx1*x1i + psx0*(x1i+x2i)*0.5+ psyy*yyi;
			pdni = pdx1*x1i + pdx0*(x1i+x2i)*0.5+ pdyy*yyi;
			psi_ni = psi_ni + qopi*(psni*ssum + pdni*sdif);
			
			p.qtanm = p.qtanm + qopi*(psni*ssum + pdni*sdif);
			
			p.dqdm(jm) = p.dqdm(jm) + qopi*(-psni*dsim + pdni*dsim);
			p.dqdm(jo) = p.dqdm(jo) + qopi*(-psni*dsio - pdni*dsio);
			p.dqdm(jp) = p.dqdm(jp) + qopi*( psni*(dsio+dsim)+ pdni*(dsio-dsim));
			
			%//------- calculate source contribution to psi	for  0-2  half-panel
			dxinv = 1.0/(x0-x2);
			psum = x2*(t2-apan) - x0*(t0-apan) + 0.5*yy*(g0-g2);
			pdif = ((x0+x2)*psum + rs0*(t0-apan) - rs2*(t2-apan)				+ (x2-x0)*yy) * dxinv;
			
			psx0 =  -(t0-apan);
			psx2 =	t2-apan;
			psyy =  0.5*(g0-g2);
			
			pdx0 = ((x0+x2)*psx0 + psum + 2.0*x0*(t0-apan) - pdif) * dxinv;
			pdx2 = ((x0+x2)*psx2 + psum - 2.0*x2*(t2-apan) + pdif) * dxinv;
			pdyy = ((x0+x2)*psyy + 2.0*(x2-x0 + yy*(t0-t2))	  ) * dxinv;
			
			dsp = sqrt((p.x(jq)-p.x(jo))*(p.x(jq)-p.x(jo)) + (p.y(jq)-p.y(jo))*(p.y(jq)-p.y(jo)));
			dsip = 1.0/dsp;
			
			%//%//ccc 	  sig2 = (p.sig(jq) - p.sig(jo))*dsip
			%//%//ccc 	  sig0 = (p.sig(jp) - p.sig(jo))*dsio
			%//%//ccc 	  ssum = sig2 + sig0
			%//%//ccc 	  sdif = sig2 - sig0
			
			ssum = (p.sig(jq) - p.sig(jo))*dsip + (p.sig(jp) - p.sig(jo))*dsio;
			sdif = (p.sig(jq) - p.sig(jo))*dsip - (p.sig(jp) - p.sig(jo))*dsio;
			
			psi = psi + qopi*(psum*ssum + pdif*sdif);
			
			%//------- dpsi/dm
			p.dzdm(jo) = p.dzdm(jo) + qopi*(-psum*(dsip+dsio)- pdif*(dsip-dsio));
			p.dzdm(jp) = p.dzdm(jp) + qopi*( psum*dsio - pdif*dsio);
			p.dzdm(jq) = p.dzdm(jq) + qopi*( psum*dsip + pdif*dsip);
			
			%//------- dpsi/dni
			psni = psx0*(x1i+x2i)*0.5 + psx2*x2i + psyy*yyi;
			pdni = pdx0*(x1i+x2i)*0.5 + pdx2*x2i + pdyy*yyi;
			psi_ni = psi_ni + qopi*(psni*ssum + pdni*sdif);
			
			p.qtanm = p.qtanm + qopi*(psni*ssum + pdni*sdif);
			
			p.dqdm(jo) = p.dqdm(jo) + qopi*(-psni*(dsip+dsio)- pdni*(dsip-dsio));
			p.dqdm(jp) = p.dqdm(jp) + qopi*( psni*dsio - pdni*dsio);
			p.dqdm(jq) = p.dqdm(jq) + qopi*( psni*dsip + pdni*dsip);
			
        end
		
		%//------ calculate vortex panel contribution to psi
		dxinv = 1.0/(x1-x2);
		
		psis = 0.5*x1*g1 - 0.5*x2*g2 + x2 - x1 + yy*(t1-t2);
		
		psid = ((x1+x2)*psis + 0.5*(rs2*g2-rs1*g1 + x1*x1-x2*x2))*dxinv;
		
		psx1 = 0.5*g1;
		psx2 = -.5*g2;
		psyy = t1-t2;
		
		pdx1 = ((x1+x2)*psx1 + psis - x1*g1 - psid)*dxinv;
		pdx2 = ((x1+x2)*psx2 + psis + x2*g2 + psid)*dxinv;
		pdyy = ((x1+x2)*psyy - yy*(g1-g2) 	  )*dxinv;
		
		gsum1 = p.gamu(jp,1) + p.gamu(jo,1);
		gsum2 = p.gamu(jp,2) + p.gamu(jo,2);
		gdif1 = p.gamu(jp,1) - p.gamu(jo,1);
		gdif2 = p.gamu(jp,2) - p.gamu(jo,2);
		
		gsum = p.gam(jp) + p.gam(jo);
		gdif = p.gam(jp) - p.gam(jo);
		
		psi = psi - qopi*(psis*gsum + psid*gdif);
		
		%//------ dpsi/dgam
		p.dzdg(jo) = p.dzdg(jo) - qopi*(psis-psid);
		p.dzdg(jp) = p.dzdg(jp) - qopi*(psis+psid);
		
		%//------ dpsi/dni
		psni = psx1*x1i + psx2*x2i + psyy*yyi;
		pdni = pdx1*x1i + pdx2*x2i + pdyy*yyi;
		psi_ni = psi_ni - qopi*(gsum*psni + gdif*pdni);
		
		p.qtan1 = p.qtan1 - qopi*(gsum1*psni + gdif1*pdni);
		p.qtan2 = p.qtan2 - qopi*(gsum2*psni + gdif2*pdni);
		
		p.dqdg(jo) = p.dqdg(jo) - qopi*(psni - pdni);
		p.dqdg(jp) = p.dqdg(jp) - qopi*(psni + pdni);
		
		if(geolin) 
			
			%//------- dpsi/dn
			p.dzdn(jo) = p.dzdn(jo)- qopi*gsum*(psx1*x1o + psx2*x2o + psyy*yyo)...
				- qopi*gdif*(pdx1*x1o + pdx2*x2o + pdyy*yyo);
			p.dzdn(jp) = p.dzdn(jp)- qopi*gsum*(psx1*x1p + psx2*x2p + psyy*yyp)...
				- qopi*gdif*(pdx1*x1p + pdx2*x2p + pdyy*yyp);
			%//------- dpsi/dp
			p.z_qdof0 = p.z_qdof0- qopi*((psis-psid)*qf0(jo) + (psis+psid)*qf0(jp));
			p.z_qdof1 = p.z_qdof1- qopi*((psis-psid)*qf1(jo) + (psis+psid)*qf1(jp));
			p.z_qdof2 = p.z_qdof2- qopi*((psis-psid)*qf2(jo) + (psis+psid)*qf2(jp));
			p.z_qdof3 = p.z_qdof3- qopi*((psis-psid)*qf3(jo) + (psis+psid)*qf3(jp));
        end
    end
%stop21:
	psig = 0.5*yy*(g1-g2) + x2*(t2-apan) - x1*(t1-apan);
	pgam = 0.5*x1*g1 - 0.5*x2*g2 + x2 - x1 + yy*(t1-t2);
	
	psigx1 = -(t1-apan);
	psigx2 =	 t2-apan;
	psigyy = 0.5*(g1-g2);
	pgamx1 = 0.5*g1;
	pgamx2 = -.5*g2;
	pgamyy = t1-t2;
	
	psigni = psigx1*x1i + psigx2*x2i + psigyy*yyi;
	pgamni = pgamx1*x1i + pgamx2*x2i + pgamyy*yyi;
	
	%//---- te panel source and vortex strengths
	sigte1 = 0.5*scs*(p.gamu(jp,1) - p.gamu(jo,1));
	sigte2 = 0.5*scs*(p.gamu(jp,2) - p.gamu(jo,2));
	gamte1 = -.5*sds*(p.gamu(jp,1) - p.gamu(jo,1));
	gamte2 = -.5*sds*(p.gamu(jp,2) - p.gamu(jo,2));
	
	sigte = 0.5*scs*(p.gam(jp) - p.gam(jo));
	gamte = -.5*sds*(p.gam(jp) - p.gam(jo));
	
	%//---- te panel contribution to psi
	psi = psi + hopi*(psig*sigte - pgam*gamte);
	
	%//---- dpsi/dgam
	p.dzdg(jo) = p.dzdg(jo) - hopi*psig*scs*0.5;
	p.dzdg(jp) = p.dzdg(jp) + hopi*psig*scs*0.5;
	
	p.dzdg(jo) = p.dzdg(jo) - hopi*pgam*sds*0.5;
	p.dzdg(jp) = p.dzdg(jp) + hopi*pgam*sds*0.5;
	
	%//---- dpsi/dni
	psi_ni = psi_ni + hopi*(psigni*sigte - pgamni*gamte);
	
	p.qtan1 = p.qtan1 + hopi*(psigni*sigte1 - pgamni*gamte1);
	p.qtan2 = p.qtan2 + hopi*(psigni*sigte2 - pgamni*gamte2);
	
	p.dqdg(jo) = p.dqdg(jo) - hopi*(psigni*0.5*scs + pgamni*0.5*sds);
	p.dqdg(jp) = p.dqdg(jp) + hopi*(psigni*0.5*scs + pgamni*0.5*sds);
	
	if(geolin) 
		
		%//----- dpsi/dn
		p.dzdn(jo) = p.dzdn(jo)			+ hopi*(psigx1*x1o + psigx2*x2o + psigyy*yyo)*sigte			- hopi*(pgamx1*x1o + pgamx2*x2o + pgamyy*yyo)*gamte;
		p.dzdn(jp) = p.dzdn(jp)			+ hopi*(psigx1*x1p + psigx2*x2p + psigyy*yyp)*sigte			- hopi*(pgamx1*x1p + pgamx2*x2p + pgamyy*yyp)*gamte;
		
		%//----- dpsi/dp
		p.z_qdof0 = p.z_qdof0 + hopi*psig*0.5*(qf0(jp)-qf0(jo))*scs			+ hopi*pgam*0.5*(qf0(jp)-qf0(jo))*sds;
		p.z_qdof1 = p.z_qdof1 + hopi*psig*0.5*(qf1(jp)-qf1(jo))*scs			+ hopi*pgam*0.5*(qf1(jp)-qf1(jo))*sds;
		p.z_qdof2 = p.z_qdof2 + hopi*psig*0.5*(qf2(jp)-qf2(jo))*scs			+ hopi*pgam*0.5*(qf2(jp)-qf2(jo))*sds;
		p.z_qdof3 = p.z_qdof3 + hopi*psig*0.5*(qf3(jp)-qf3(jo))*scs			+ hopi*pgam*0.5*(qf3(jp)-qf3(jo))*sds;
		

    end
return
                   
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    end
                end
            end
        end
		
		dso = sqrt((p.x(jo)-p.x(jp))*(p.x(jo)-p.x(jp)) + (p.y(jo)-p.y(jp))*(p.y(jo)-p.y(jp)));
		
		%//------ skip null panel
		if(dso == 0.0) 
            %goto stop10; %//check - unsafe comparison
		    nothing=1;
            continue;
        end
		dsio = 1.0 / dso;
		
		apan = p.apanel(jo);
		
		rx1 = xi - p.x(jo);
		ry1 = yi - p.y(jo);
		rx2 = xi - p.x(jp);
		ry2 = yi - p.y(jp);
		
		sx = (p.x(jp) - p.x(jo)) * dsio;
		sy = (p.y(jp) - p.y(jo)) * dsio;
		
		x1 = sx*rx1 + sy*ry1;
		x2 = sx*rx2 + sy*ry2;
		yy = sx*ry1 - sy*rx1;
		
		rs1 = rx1*rx1 + ry1*ry1;
		rs2 = rx2*rx2 + ry2*ry2;
		
		%//------ set reflection flag sgn to avoid branch problems with arctan
		if(io>=1 && io<=p.n) 
			%//------- no problem on airfoil surface
			sgn = 1.0;
		
        else
			%//------- make sure arctan falls between  -/+  pi/2
			sgn = 1*sign(yy);
        end
		
		%//------ set log(r^2) and arctan(p.x/p.y), correcting for reflection if any
		if(io~=jo && rs1>0.0) 
			g1 = log(rs1);
			t1 = atan2(sgn*x1,sgn*yy) + (0.5- 0.5*sgn)*pi;
		
        else
			g1 = 0.0;
			t1 = 0.0;
        end
		
		if(io~=jp && rs2>0.0) 
			g2 = log(rs2);
			t2 = atan2(sgn*x2,sgn*yy) + (0.5- 0.5*sgn)*pi;
		
        else
			g2 = 0.0;
			t2 = 0.0;
        end
		
		x1i = sx*nxi + sy*nyi;
		x2i = sx*nxi + sy*nyi;
		yyi = sx*nyi - sy*nxi;
		
		if(geolin) 
			nxo = p.nx(jo);
			nyo = p.ny(jo);
			nxp = p.nx(jp);
			nyp = p.ny(jp);
			
			x1o =-((rx1-x1*sx)*nxo + (ry1-x1*sy)*nyo)*dsio-(sx*nxo+sy*nyo);
			x1p = ((rx1-x1*sx)*nxp + (ry1-x1*sy)*nyp)*dsio;
			x2o =-((rx2-x2*sx)*nxo + (ry2-x2*sy)*nyo)*dsio;
			x2p = ((rx2-x2*sx)*nxp + (ry2-x2*sy)*nyp)*dsio-(sx*nxp+sy*nyp);
			yyo = ((rx1+x1*sy)*nyo - (ry1-x1*sx)*nxo)*dsio-(sx*nyo-sy*nxo);
			yyp =-((rx1-x1*sy)*nyp - (ry1+x1*sx)*nxp)*dsio;
        end
		
		if (jo==p.n) 
            continue;
        end
		
		if(siglin) 
			
			%//------- set up midpo quantities
			x0 = 0.5*(x1+x2);
			rs0 = x0*x0 + yy*yy;
			g0 = log(rs0);
			t0 = atan2(sgn*x0,sgn*yy) + (0.5- 0.5*sgn)*pi;
			
			%//------- calculate source contribution to psi	for  1-0  half-panel
			dxinv = 1.0/(x1-x0);
			psum = x0*(t0-apan) - x1*(t1-apan) + 0.5*yy*(g1-g0);
			pdif = ((x1+x0)*psum + rs1*(t1-apan) - rs0*(t0-apan)+ (x0-x1)*yy) * dxinv;
			
			psx1 =	-(t1-apan);
			psx0 =	  t0-apan;
			psyy =	0.5*(g1-g0);
			
			pdx1 = ((x1+x0)*psx1 + psum + 2.0*x1*(t1-apan) - pdif) * dxinv;
			pdx0 = ((x1+x0)*psx0 + psum - 2.0*x0*(t0-apan) + pdif) * dxinv;
			pdyy = ((x1+x0)*psyy + 2.0*(x0-x1 + yy*(t1-t0))   ) * dxinv;
			
			dsm = sqrt((p.x(jp)-p.x(jm))*(p.x(jp)-p.x(jm)) + (p.y(jp)-p.y(jm))*(p.y(jp)-p.y(jm)));
			dsim = 1.0/dsm;
			
			%//%//ccc	   sig0 = (p.sig(jp) - p.sig(jo))*dsio
			%//%//ccc	   sig1 = (p.sig(jp) - p.sig(jm))*dsim
			%//%//ccc	   ssum = sig0 + sig1
			%//%//ccc	   sdif = sig0 - sig1 
			
			ssum = (p.sig(jp) - p.sig(jo))*dsio + (p.sig(jp) - p.sig(jm))*dsim;
			sdif = (p.sig(jp) - p.sig(jo))*dsio - (p.sig(jp) - p.sig(jm))*dsim;
			
			psi = psi + qopi*(psum*ssum + pdif*sdif);
			
			%//------- dpsi/dm
			p.dzdm(jm) = p.dzdm(jm) + qopi*(-psum*dsim + pdif*dsim);
			p.dzdm(jo) = p.dzdm(jo) + qopi*(-psum*dsio - pdif*dsio);
			p.dzdm(jp) = p.dzdm(jp) + qopi*( psum*(dsio+dsim) + pdif*(dsio-dsim));
			
			%//------- dpsi/dni
			psni = psx1*x1i + psx0*(x1i+x2i)*0.5 + psyy*yyi;
			pdni = pdx1*x1i + pdx0*(x1i+x2i)*0.5 + pdyy*yyi;
			psi_ni = psi_ni + qopi*(psni*ssum + pdni*sdif);
			
			p.qtanm = p.qtanm + qopi*(psni*ssum + pdni*sdif);
			
			p.dqdm(jm) = p.dqdm(jm) + qopi*(-psni*dsim + pdni*dsim);
			p.dqdm(jo) = p.dqdm(jo) + qopi*(-psni*dsio - pdni*dsio);
			p.dqdm(jp) = p.dqdm(jp) + qopi*( psni*(dsio+dsim)+ pdni*(dsio-dsim));
			
			
			%//------- calculate source contribution to psi	for  0-2  half-panel
			dxinv = 1.0/(x0-x2);
			psum = x2*(t2-apan) - x0*(t0-apan) + 0.5*yy*(g0-g2);
			pdif = ((x0+x2)*psum + rs0*(t0-apan) - rs2*(t2-apan)+ (x2-x0)*yy) * dxinv;
			
			psx0 =  -(t0-apan);
			psx2 =	  t2-apan;
			psyy =	0.5*(g0-g2);
			
			pdx0 = ((x0+x2)*psx0 + psum + 2.0*x0*(t0-apan) - pdif) * dxinv;
			pdx2 = ((x0+x2)*psx2 + psum - 2.0*x2*(t2-apan) + pdif) * dxinv;
			pdyy = ((x0+x2)*psyy + 2.0*(x2-x0 + yy*(t0-t2))	  ) * dxinv;
			
			dsp = sqrt((p.x(jq)-p.x(jo))*(p.x(jq)-p.x(jo)) + (p.y(jq)-p.y(jo))*(p.y(jq)-p.y(jo)));
			dsip = 1.0/dsp;
			
			%//%//ccc		  sig2 = (p.sig(jq) - p.sig(jo))*dsip
			%//%//ccc		  sig0 = (p.sig(jp) - p.sig(jo))*dsio
			%//%//ccc		  ssum = sig2 + sig0
			%//%//ccc		  sdif = sig2 - sig0
			
			ssum = (p.sig(jq) - p.sig(jo))*dsip + (p.sig(jp) - p.sig(jo))*dsio;
			sdif = (p.sig(jq) - p.sig(jo))*dsip - (p.sig(jp) - p.sig(jo))*dsio;
			
			psi = psi + qopi*(psum*ssum + pdif*sdif);
			
			%//------- dpsi/dm
			p.dzdm(jo) = p.dzdm(jo) + qopi*(-psum*(dsip+dsio)- pdif*(dsip-dsio));
			p.dzdm(jp) = p.dzdm(jp) + qopi*( psum*dsio - pdif*dsio);
			p.dzdm(jq) = p.dzdm(jq) + qopi*( psum*dsip + pdif*dsip);
			
			%//------- dpsi/dni
			psni = psx0*(x1i+x2i)*0.5 + psx2*x2i + psyy*yyi;
			pdni = pdx0*(x1i+x2i)*0.5 + pdx2*x2i + pdyy*yyi;
			psi_ni = psi_ni + qopi*(psni*ssum + pdni*sdif);
			
			p.qtanm = p.qtanm + qopi*(psni*ssum + pdni*sdif);
			
			p.dqdm(jo) = p.dqdm(jo) + qopi*(-psni*(dsip+dsio)- pdni*(dsip-dsio));
			p.dqdm(jp) = p.dqdm(jp) + qopi*( psni*dsio - pdni*dsio);
			p.dqdm(jq) = p.dqdm(jq) + qopi*( psni*dsip + pdni*dsip);
			
        end
		
		%//------ calculate vortex panel contribution to psi
		dxinv = 1.0/(x1-x2);
		psis = 0.5*x1*g1 - 0.5*x2*g2 + x2 - x1 + yy*(t1-t2);
		psid = ((x1+x2)*psis + 0.5*(rs2*g2-rs1*g1 + x1*x1-x2*x2))*dxinv;
		
		psx1 = 0.5*g1;
		psx2 = -.5*g2;
		psyy = t1-t2;
		
		pdx1 = ((x1+x2)*psx1 + psis - x1*g1 - psid)*dxinv;
		pdx2 = ((x1+x2)*psx2 + psis + x2*g2 + psid)*dxinv;
		pdyy = ((x1+x2)*psyy - yy*(g1-g2)		  )*dxinv;
		
		gsum1 = p.gamu(jp,1) + p.gamu(jo,1);
		gsum2 = p.gamu(jp,2) + p.gamu(jo,2);
		gdif1 = p.gamu(jp,1) - p.gamu(jo,1);
		gdif2 = p.gamu(jp,2) - p.gamu(jo,2);
		
		gsum = p.gam(jp) + p.gam(jo);
		gdif = p.gam(jp) - p.gam(jo);
		
		psi = psi + qopi*(psis*gsum + psid*gdif);
		
		%//------ dpsi/dgam
		p.dzdg(jo) = p.dzdg(jo) + qopi*(psis-psid);
		p.dzdg(jp) = p.dzdg(jp) + qopi*(psis+psid);
		
		%//------ dpsi/dni
		psni = psx1*x1i + psx2*x2i + psyy*yyi;
		pdni = pdx1*x1i + pdx2*x2i + pdyy*yyi;
		psi_ni = psi_ni + qopi*(gsum*psni + gdif*pdni);
		
		p.qtan1 = p.qtan1 + qopi*(gsum1*psni + gdif1*pdni);
		p.qtan2 = p.qtan2 + qopi*(gsum2*psni + gdif2*pdni);
		
		p.dqdg(jo) = p.dqdg(jo) + qopi*(psni - pdni);
		p.dqdg(jp) = p.dqdg(jp) + qopi*(psni + pdni);
		
		if(geolin) 
			
			%//------- dpsi/dn
			p.dzdn(jo) = p.dzdn(jo)+ qopi*gsum*(psx1*x1o + psx2*x2o + psyy*yyo)				+ qopi*gdif*(pdx1*x1o + pdx2*x2o + pdyy*yyo);
			p.dzdn(jp) = p.dzdn(jp)+ qopi*gsum*(psx1*x1p + psx2*x2p + psyy*yyp)				+ qopi*gdif*(pdx1*x1p + pdx2*x2p + pdyy*yyp);
			%//------- dpsi/dp
			p.z_qdof0 = p.z_qdof0				+ qopi*((psis-psid)*qf0(jo) + (psis+psid)*qf0(jp));
			p.z_qdof1 = p.z_qdof1				+ qopi*((psis-psid)*qf1(jo) + (psis+psid)*qf1(jp));
			p.z_qdof2 = p.z_qdof2				+ qopi*((psis-psid)*qf2(jo) + (psis+psid)*qf2(jp));
			p.z_qdof3 = p.z_qdof3				+ qopi*((psis-psid)*qf3(jo) + (psis+psid)*qf3(jp));
        end
%stop10:
		nothing = 1;
    end

%stop11:
	psig = 0.5*yy*(g1-g2) + x2*(t2-apan) - x1*(t1-apan);
	pgam = 0.5*x1*g1 - 0.5*x2*g2 + x2 - x1 + yy*(t1-t2);
	
	psigx1 = -(t1-apan);
	psigx2 =	 t2-apan;
	psigyy = 0.5*(g1-g2);
	pgamx1 = 0.5*g1;
	pgamx2 = -.5*g2;
	pgamyy = t1-t2;
	
	psigni = psigx1*x1i + psigx2*x2i + psigyy*yyi;
	pgamni = pgamx1*x1i + pgamx2*x2i + pgamyy*yyi;
	
	%//---- te panel source and vortex strengths
	sigte1 = 0.5*scs*(p.gamu(jp,1) - p.gamu(jo,1));
	sigte2 = 0.5*scs*(p.gamu(jp,2) - p.gamu(jo,2));
	gamte1 = -.5*sds*(p.gamu(jp,1) - p.gamu(jo,1));
	gamte2 = -.5*sds*(p.gamu(jp,2) - p.gamu(jo,2));
	
	sigte = 0.5*scs*(p.gam(jp) - p.gam(jo));
	gamte = -.5*sds*(p.gam(jp) - p.gam(jo));
	
	%//---- te panel contribution to psi
	psi = psi + hopi*(psig*sigte + pgam*gamte);
	
	%//---- dpsi/dgam
	p.dzdg(jo) = p.dzdg(jo) - hopi*psig*scs*0.5;
	p.dzdg(jp) = p.dzdg(jp) + hopi*psig*scs*0.5;
	
	p.dzdg(jo) = p.dzdg(jo) + hopi*pgam*sds*0.5;
	p.dzdg(jp) = p.dzdg(jp) - hopi*pgam*sds*0.5;
	
	%//---- dpsi/dni
	psi_ni = psi_ni + hopi*(psigni*sigte + pgamni*gamte);
	
	p.qtan1 = p.qtan1 + hopi*(psigni*sigte1 + pgamni*gamte1);
	p.qtan2 = p.qtan2 + hopi*(psigni*sigte2 + pgamni*gamte2);
	
	p.dqdg(jo) = p.dqdg(jo) - hopi*(psigni*0.5*scs - pgamni*0.5*sds);
	p.dqdg(jp) = p.dqdg(jp) + hopi*(psigni*0.5*scs - pgamni*0.5*sds);
	
	if(geolin)
		
		%//----- dpsi/dn
		p.dzdn(jo) = p.dzdn(jo)			+ hopi*(psigx1*x1o + psigx2*x2o + psigyy*yyo)*sigte			+ hopi*(pgamx1*x1o + pgamx2*x2o + pgamyy*yyo)*gamte;
		p.dzdn(jp) = p.dzdn(jp)			+ hopi*(psigx1*x1p + psigx2*x2p + psigyy*yyp)*sigte			+ hopi*(pgamx1*x1p + pgamx2*x2p + pgamyy*yyp)*gamte;
		
		%//----- dpsi/dp
		p.z_qdof0 = p.z_qdof0 + hopi*psig*0.5*(qf0(jp)-qf0(jo))*scs			- hopi*pgam*0.5*(qf0(jp)-qf0(jo))*sds;
		p.z_qdof1 = p.z_qdof1 + hopi*psig*0.5*(qf1(jp)-qf1(jo))*scs			- hopi*pgam*0.5*(qf1(jp)-qf1(jo))*sds;
		p.z_qdof2 = p.z_qdof2 + hopi*psig*0.5*(qf2(jp)-qf2(jo))*scs			- hopi*pgam*0.5*(qf2(jp)-qf2(jo))*sds;
		p.z_qdof3 = p.z_qdof3 + hopi*psig*0.5*(qf3(jp)-qf3(jo))*scs			- hopi*pgam*0.5*(qf3(jp)-qf3(jo))*sds;
		
    end   
                        	%//**** freestream terms
	psi = psi + qinf*(cosa*yi - sina*xi);
	
	%//---- dpsi/dn
	psi_ni = psi_ni + qinf*(cosa*nyi - sina*nxi);
	
	p.qtan1 = p.qtan1 + qinf*nyi;
	p.qtan2 = p.qtan2 - qinf*nxi;
	
	%//---- dpsi/dqinf
	p.z_qinf = p.z_qinf + (cosa*yi - sina*xi);
	
	%//---- dpsi/dalfa
	p.z_alfa = p.z_alfa - qinf*(sina*yi + cosa*xi);
	
	if(limage==false)
		return;
    end

	for(jo=1:p.n)
		jp = jo+1;
		
		jm = jo-1;
		jq = jp+1;
		
		if(jo==1) 
            jm = jo;
        else
			if(jo==p.n-1)  jq = jp;
			else 
				if(jo==p.n) jp = 1
				if((p.x(jo)-p.x(jp))*(p.x(jo)-p.x(jp)) + (p.y(jo)-p.y(jp))*(p.y(jo)-p.y(jp)) < seps*seps)
					return;
                end
                end
            end
        end
		
		dso = sqrt((p.x(jo)-p.x(jp))*(p.x(jo)-p.x(jp)) + (p.y(jo)-p.y(jp))*(p.y(jo)-p.y(jp)));
		
		%//------ skip null panel
		if (dso == 0.0) 
            %goto stop20; %// check - unsafe condition
            
            		nothing = 1;
                    continue;
        end

            
    
            
            
            
            
		dsio = 1.0 / dso;
		%//%//ccc   apan = p.apanel(jo)
		apan = pi - p.apanel(jo) + 2.0*p.alfa;
		
		xjo = p.x(jo) + 2.0*(yimage+p.y(jo))*sina;
		yjo = p.y(jo) - 2.0*(yimage+p.y(jo))*cosa;
		xjp = p.x(jp) + 2.0*(yimage+p.y(jp))*sina;
		yjp = p.y(jp) - 2.0*(yimage+p.y(jp))*cosa;
		
		rx1 = xi - xjo;
		ry1 = yi - yjo;
		rx2 = xi - xjp;
		ry2 = yi - yjp;
		
		sx = (xjp - xjo) * dsio;
		sy = (yjp - yjo) * dsio;
		
		x1 = sx*rx1 + sy*ry1;
		x2 = sx*rx2 + sy*ry2;
		yy = sx*ry1 - sy*rx1;
		
		rs1 = rx1*rx1 + ry1*ry1;
		rs2 = rx2*rx2 + ry2*ry2;
		
		%//------ set reflection flag sgn to avoid branch problems with arctan
		if(io>=1 && io<=p.n) 
			%//------- no problem on airfoil surface
			sgn = 1.0;
		
        else
			%//------- make sure arctan falls between	-/+  pi/2
			sgn = 1*sign(yy);
        end
		
		%//------ set log(r^2) and arctan(p.x/p.y), correcting for reflection if any
		g1 = log(rs1);
		t1 = atan2(sgn*x1,sgn*yy) + (0.5 - 0.5*sgn)*pi;
		
		g2 = log(rs2);
		t2 = atan2(sgn*x2,sgn*yy) + (0.5 - 0.5*sgn)*pi;
		
		x1i = sx*nxi + sy*nyi;
		x2i = sx*nxi + sy*nyi;
		yyi = sx*nyi - sy*nxi;
		
		if(geolin) 
			nxo = p.nx(jo);
			nyo = p.ny(jo);
			nxp = p.nx(jp);
			nyp = p.ny(jp);
			
			x1o =-((rx1-x1*sx)*nxo + (ry1-x1*sy)*nyo)*dsio-(sx*nxo+sy*nyo);
			x1p = ((rx1-x1*sx)*nxp + (ry1-x1*sy)*nyp)*dsio;
			x2o =-((rx2-x2*sx)*nxo + (ry2-x2*sy)*nyo)*dsio;
			x2p = ((rx2-x2*sx)*nxp + (ry2-x2*sy)*nyp)*dsio-(sx*nxp+sy*nyp);
			yyo = ((rx1+x1*sy)*nyo - (ry1-x1*sx)*nxo)*dsio-(sx*nyo-sy*nxo);
			yyp =-((rx1-x1*sy)*nyp - (ry1+x1*sx)*nxp)*dsio;
        end
		
		if(jo==p.n) 
            break;
        end
			
		if(siglin) 
			
			%//------- set up midpo quantities
			x0 = 0.5*(x1+x2);
			rs0 = x0*x0 + yy*yy;
			g0 = log(rs0);
			t0 = atan2(sgn*x0,sgn*yy) + (0.5 - 0.5*sgn)*pi;
			
			%//------- calculate source contribution to psi	for  1-0  half-panel
			dxinv = 1.0/(x1-x0);
			psum = x0*(t0-apan) - x1*(t1-apan) + 0.5*yy*(g1-g0);
			pdif = ((x1+x0)*psum + rs1*(t1-apan) - rs0*(t0-apan)				+ (x0-x1)*yy) * dxinv;
			
			psx1 =  -(t1-apan);
			psx0 =	t0-apan;
			psyy =  0.5*(g1-g0);
			
			pdx1 = ((x1+x0)*psx1 + psum + 2.0*x1*(t1-apan) - pdif) * dxinv;
			pdx0 = ((x1+x0)*psx0 + psum - 2.0*x0*(t0-apan) + pdif) * dxinv;
			pdyy = ((x1+x0)*psyy + 2.0*(x0-x1 + yy*(t1-t0))	  ) * dxinv;
			
			dsm = sqrt((p.x(jp)-p.x(jm))*(p.x(jp)-p.x(jm)) + (p.y(jp)-p.y(jm))*(p.y(jp)-p.y(jm)));
			dsim = 1.0/dsm;
			
			%//%//ccc    sig0 = (p.sig(jp) - p.sig(jo))*dsio
			%//%//ccc    sig1 = (p.sig(jp) - p.sig(jm))*dsim
			%//%//ccc    ssum = sig0 + sig1
			%//%//ccc    sdif = sig0 - sig1
			
			ssum = (p.sig(jp) - p.sig(jo))*dsio + (p.sig(jp) - p.sig(jm))*dsim;
			sdif = (p.sig(jp) - p.sig(jo))*dsio - (p.sig(jp) - p.sig(jm))*dsim;
			
			psi = psi + qopi*(psum*ssum + pdif*sdif);
			
			%//------- dpsi/dm
			p.dzdm(jm) = p.dzdm(jm) + qopi*(-psum*dsim + pdif*dsim);
			p.dzdm(jo) = p.dzdm(jo) + qopi*(-psum*dsio - pdif*dsio);
			p.dzdm(jp) = p.dzdm(jp) + qopi*( psum*(dsio+dsim)+ pdif*(dsio-dsim));
			
			%//------- dpsi/dni
			psni = psx1*x1i + psx0*(x1i+x2i)*0.5+ psyy*yyi;
			pdni = pdx1*x1i + pdx0*(x1i+x2i)*0.5+ pdyy*yyi;
			psi_ni = psi_ni + qopi*(psni*ssum + pdni*sdif);
			
			p.qtanm = p.qtanm + qopi*(psni*ssum + pdni*sdif);
			
			p.dqdm(jm) = p.dqdm(jm) + qopi*(-psni*dsim + pdni*dsim);
			p.dqdm(jo) = p.dqdm(jo) + qopi*(-psni*dsio - pdni*dsio);
			p.dqdm(jp) = p.dqdm(jp) + qopi*( psni*(dsio+dsim)+ pdni*(dsio-dsim));
			
			%//------- calculate source contribution to psi	for  0-2  half-panel
			dxinv = 1.0/(x0-x2);
			psum = x2*(t2-apan) - x0*(t0-apan) + 0.5*yy*(g0-g2);
			pdif = ((x0+x2)*psum + rs0*(t0-apan) - rs2*(t2-apan)				+ (x2-x0)*yy) * dxinv;
			
			psx0 =  -(t0-apan);
			psx2 =	t2-apan;
			psyy =  0.5*(g0-g2);
			
			pdx0 = ((x0+x2)*psx0 + psum + 2.0*x0*(t0-apan) - pdif) * dxinv;
			pdx2 = ((x0+x2)*psx2 + psum - 2.0*x2*(t2-apan) + pdif) * dxinv;
			pdyy = ((x0+x2)*psyy + 2.0*(x2-x0 + yy*(t0-t2))	  ) * dxinv;
			
			dsp = sqrt((p.x(jq)-p.x(jo))*(p.x(jq)-p.x(jo)) + (p.y(jq)-p.y(jo))*(p.y(jq)-p.y(jo)));
			dsip = 1.0/dsp;
			
			%//%//ccc 	  sig2 = (p.sig(jq) - p.sig(jo))*dsip
			%//%//ccc 	  sig0 = (p.sig(jp) - p.sig(jo))*dsio
			%//%//ccc 	  ssum = sig2 + sig0
			%//%//ccc 	  sdif = sig2 - sig0
			
			ssum = (p.sig(jq) - p.sig(jo))*dsip + (p.sig(jp) - p.sig(jo))*dsio;
			sdif = (p.sig(jq) - p.sig(jo))*dsip - (p.sig(jp) - p.sig(jo))*dsio;
			
			psi = psi + qopi*(psum*ssum + pdif*sdif);
			
			%//------- dpsi/dm
			p.dzdm(jo) = p.dzdm(jo) + qopi*(-psum*(dsip+dsio)- pdif*(dsip-dsio));
			p.dzdm(jp) = p.dzdm(jp) + qopi*( psum*dsio - pdif*dsio);
			p.dzdm(jq) = p.dzdm(jq) + qopi*( psum*dsip + pdif*dsip);
			
			%//------- dpsi/dni
			psni = psx0*(x1i+x2i)*0.5 + psx2*x2i + psyy*yyi;
			pdni = pdx0*(x1i+x2i)*0.5 + pdx2*x2i + pdyy*yyi;
			psi_ni = psi_ni + qopi*(psni*ssum + pdni*sdif);
			
			p.qtanm = p.qtanm + qopi*(psni*ssum + pdni*sdif);
			
			p.dqdm(jo) = p.dqdm(jo) + qopi*(-psni*(dsip+dsio)- pdni*(dsip-dsio));
			p.dqdm(jp) = p.dqdm(jp) + qopi*( psni*dsio - pdni*dsio);
			p.dqdm(jq) = p.dqdm(jq) + qopi*( psni*dsip + pdni*dsip);
			
        end
		
		%//------ calculate vortex panel contribution to psi
		dxinv = 1.0/(x1-x2);
		
		psis = 0.5*x1*g1 - 0.5*x2*g2 + x2 - x1 + yy*(t1-t2);
		
		psid = ((x1+x2)*psis + 0.5*(rs2*g2-rs1*g1 + x1*x1-x2*x2))*dxinv;
		
		psx1 = 0.5*g1;
		psx2 = -.5*g2;
		psyy = t1-t2;
		
		pdx1 = ((x1+x2)*psx1 + psis - x1*g1 - psid)*dxinv;
		pdx2 = ((x1+x2)*psx2 + psis + x2*g2 + psid)*dxinv;
		pdyy = ((x1+x2)*psyy - yy*(g1-g2) 	  )*dxinv;
		
		gsum1 = p.gamu(jp,1) + p.gamu(jo,1);
		gsum2 = p.gamu(jp,2) + p.gamu(jo,2);
		gdif1 = p.gamu(jp,1) - p.gamu(jo,1);
		gdif2 = p.gamu(jp,2) - p.gamu(jo,2);
		
		gsum = p.gam(jp) + p.gam(jo);
		gdif = p.gam(jp) - p.gam(jo);
		
		psi = psi - qopi*(psis*gsum + psid*gdif);
		
		%//------ dpsi/dgam
		p.dzdg(jo) = p.dzdg(jo) - qopi*(psis-psid);
		p.dzdg(jp) = p.dzdg(jp) - qopi*(psis+psid);
		
		%//------ dpsi/dni
		psni = psx1*x1i + psx2*x2i + psyy*yyi;
		pdni = pdx1*x1i + pdx2*x2i + pdyy*yyi;
		psi_ni = psi_ni - qopi*(gsum*psni + gdif*pdni);
		
		p.qtan1 = p.qtan1 - qopi*(gsum1*psni + gdif1*pdni);
		p.qtan2 = p.qtan2 - qopi*(gsum2*psni + gdif2*pdni);
		
		p.dqdg(jo) = p.dqdg(jo) - qopi*(psni - pdni);
		p.dqdg(jp) = p.dqdg(jp) - qopi*(psni + pdni);
		
		if(geolin) 
			
			%//------- dpsi/dn
			p.dzdn(jo) = p.dzdn(jo)- qopi*gsum*(psx1*x1o + psx2*x2o + psyy*yyo)...
				- qopi*gdif*(pdx1*x1o + pdx2*x2o + pdyy*yyo);
			p.dzdn(jp) = p.dzdn(jp)- qopi*gsum*(psx1*x1p + psx2*x2p + psyy*yyp)...
				- qopi*gdif*(pdx1*x1p + pdx2*x2p + pdyy*yyp);
			%//------- dpsi/dp
			p.z_qdof0 = p.z_qdof0- qopi*((psis-psid)*qf0(jo) + (psis+psid)*qf0(jp));
			p.z_qdof1 = p.z_qdof1- qopi*((psis-psid)*qf1(jo) + (psis+psid)*qf1(jp));
			p.z_qdof2 = p.z_qdof2- qopi*((psis-psid)*qf2(jo) + (psis+psid)*qf2(jp));
			p.z_qdof3 = p.z_qdof3- qopi*((psis-psid)*qf3(jo) + (psis+psid)*qf3(jp));
        end
    end
%stop21:
	psig = 0.5*yy*(g1-g2) + x2*(t2-apan) - x1*(t1-apan);
	pgam = 0.5*x1*g1 - 0.5*x2*g2 + x2 - x1 + yy*(t1-t2);
	
	psigx1 = -(t1-apan);
	psigx2 =	 t2-apan;
	psigyy = 0.5*(g1-g2);
	pgamx1 = 0.5*g1;
	pgamx2 = -.5*g2;
	pgamyy = t1-t2;
	
	psigni = psigx1*x1i + psigx2*x2i + psigyy*yyi;
	pgamni = pgamx1*x1i + pgamx2*x2i + pgamyy*yyi;
	
	%//---- te panel source and vortex strengths
	sigte1 = 0.5*scs*(p.gamu(jp,1) - p.gamu(jo,1));
	sigte2 = 0.5*scs*(p.gamu(jp,2) - p.gamu(jo,2));
	gamte1 = -.5*sds*(p.gamu(jp,1) - p.gamu(jo,1));
	gamte2 = -.5*sds*(p.gamu(jp,2) - p.gamu(jo,2));
	
	sigte = 0.5*scs*(p.gam(jp) - p.gam(jo));
	gamte = -.5*sds*(p.gam(jp) - p.gam(jo));
	
	%//---- te panel contribution to psi
	psi = psi + hopi*(psig*sigte - pgam*gamte);
	
	%//---- dpsi/dgam
	p.dzdg(jo) = p.dzdg(jo) - hopi*psig*scs*0.5;
	p.dzdg(jp) = p.dzdg(jp) + hopi*psig*scs*0.5;
	
	p.dzdg(jo) = p.dzdg(jo) - hopi*pgam*sds*0.5;
	p.dzdg(jp) = p.dzdg(jp) + hopi*pgam*sds*0.5;
	
	%//---- dpsi/dni
	psi_ni = psi_ni + hopi*(psigni*sigte - pgamni*gamte);
	
	p.qtan1 = p.qtan1 + hopi*(psigni*sigte1 - pgamni*gamte1);
	p.qtan2 = p.qtan2 + hopi*(psigni*sigte2 - pgamni*gamte2);
	
	p.dqdg(jo) = p.dqdg(jo) - hopi*(psigni*0.5*scs + pgamni*0.5*sds);
	p.dqdg(jp) = p.dqdg(jp) + hopi*(psigni*0.5*scs + pgamni*0.5*sds);
	
	if(geolin) 
		
		%//----- dpsi/dn
		p.dzdn(jo) = p.dzdn(jo)			+ hopi*(psigx1*x1o + psigx2*x2o + psigyy*yyo)*sigte			- hopi*(pgamx1*x1o + pgamx2*x2o + pgamyy*yyo)*gamte;
		p.dzdn(jp) = p.dzdn(jp)			+ hopi*(psigx1*x1p + psigx2*x2p + psigyy*yyp)*sigte			- hopi*(pgamx1*x1p + pgamx2*x2p + pgamyy*yyp)*gamte;
		
		%//----- dpsi/dp
		p.z_qdof0 = p.z_qdof0 + hopi*psig*0.5*(qf0(jp)-qf0(jo))*scs			+ hopi*pgam*0.5*(qf0(jp)-qf0(jo))*sds;
		p.z_qdof1 = p.z_qdof1 + hopi*psig*0.5*(qf1(jp)-qf1(jo))*scs			+ hopi*pgam*0.5*(qf1(jp)-qf1(jo))*sds;
		p.z_qdof2 = p.z_qdof2 + hopi*psig*0.5*(qf2(jp)-qf2(jo))*scs			+ hopi*pgam*0.5*(qf2(jp)-qf2(jo))*sds;
		p.z_qdof3 = p.z_qdof3 + hopi*psig*0.5*(qf3(jp)-qf3(jo))*scs			+ hopi*pgam*0.5*(qf3(jp)-qf3(jo))*sds;
		

    end
return
end
function [a, indx]=ludcmp(n, a)
% ////--------------------------------------
% //     ludcmp function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% {
% 	//    *******************************************************
% 	//    *                                                     *
% 	//    *   factors a full nxn matrix into an lu form.        *
% 	//    *   subr. baksub can back-substitute it with some rhs.*
% 	//    *   assumes matrix is non-singular...                 *
% 	//    *    ...if it isn"t, a divide by zero will result.    *
% 	//    *                                                     *
% 	//    *   a is the matrix...                                *
% 	//    *     ...replaced with its lu factors.                *
% 	//    *                                                     *
% 	//    *                              mark drela  1988       *
% 	//    *******************************************************
% 	//
% 	//
% 	bimaxok = false;
	imax =0;%//added arcds
% 	nvx=IQX;
% 	int i, j, k;
% 	double vv(IQX);
% 	double dum, sum, aamax;
% 	if(n>nvx) {
% 		CString str;
% 		str.Format("Stop ludcmp: array overflow. Increase nvx");
% 		if(m_bTrace)pXFile->WriteString(str);
% 		return false;
% 	}
	
    for (i=1:n)
		aamax = 0.0;
		for (j=1:n) 
            aamax = max(abs(a(i,j)), aamax);
        end
		vv(i) = 1.0/aamax;
    end
	
	for(j=1:n)
		for(i=1:j-1)
			sum = a(i,j);
			for (k=1:i-1) 
                sum = sum - a(i,k)*a(k,j);
            end
			a(i,j) = sum;
        end
		
		aamax = 0.0;

		for (i=j:n)
			sum = a(i,j);
			for (k=1:j-1) 
                sum = sum - a(i,k)*a(k,j);
            end
			a(i,j) = sum ;
			dum = (vv(i)*abs(sum));
			if(dum>=aamax)
				imax = i;
				aamax = dum;
				bimaxok = true;
            end
        end
% 		ASSERT(bimaxok);// to check imax has been initialized
		if(j~=imax) 
			for (k=1:n)
				dum = a(imax,k);
				a(imax,k) = a(j,k);
				a(j,k) = dum;
            end
			vv(imax) = vv(j);
        end
		
		indx(j) = imax;
		if(j~=n) 
			dum = 1.0/a(j,j);
			for(i=j+1:n)
                a(i,j) = a(i,j)*dum;
            end
        end
    end
% 	return true;
end
function [b]=baksub( n,  a,  indx,  b)
% ////--------------------------------------
% //     baksub function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% {
% 	 sum;
% 	 i, ii, ll, j;
	ii = 0;
	for (i=1:n)
		ll = indx(i);
		sum = b(ll);
		b(ll) = b(i);
		if(ii~=0) 
			for (j=ii:i-1) 
                sum = sum - a(i,j)*b(j);
            end
		else
			if(sum~=0.0) 
                ii = i;
            end
        end

		b(i) = sum;
    end
% 	//
    for (i=n:-1:1)
        sum = b(i);
        if(i<n) 
			for (j=i+1:n) 
                sum = sum - a(i,j)*b(j);
            end
        end
		
        b(i) = sum/a(i,i);
    end
% 	//
% 	return true;
end
function p=qiset(p)
% ////--------------------------------------
% //     qiset function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% 	//-------------------------------------------------------
% 	//    sets inviscid panel tangential velocity for
% 	//     current alpha.
% 	//-------------------------------------------------------
% 	
% 	p.cosa = cos(alfa);
% 	p.sina = sin(alfa);
   
	for (i=1:p.n+p.nw)
        p.qinv  (i) =  p.cosa*p.qinvu(i,1) + p.sina*p.qinvu(i,2);
        p.qinv_a(i) = -p.sina*p.qinvu(i,1) + p.cosa*p.qinvu(i,2);
    end
	
% 	return true;
	
end
function p=clcalc(xref, yref,p)
% ////--------------------------------------
% //     clcalc function
% ////--------------------------------------
% 	
%    Luca Virtuani
%    Politecnico di Milano - 2015
% //----------------------------------------

% {
% 	// modified arcds : all other variables are member variables (ex fortran common)
% 	//-----------------------------------------------------------
% 	//	   integrates surface pressures to get p.cl and p.cm.
% 	//	   integrates skin friction to get cdf.
% 	//	   calculates dcl/dalpha for prescribed-p.cl routines.
% 	//-----------------------------------------------------------
% 	
% 	//arcds addition : calculate XCp position
% 
% 	//---- moment-reference coordinates
% 	////ccc	   xref = 0.25
% 	////ccc	   yref = 0.
% 	
% 	 p.beta, p.beta_msq, p.bfac, p.bfac_msq, p.cginc;
% 	 p.cpi_gam, p.cpc_cpi;	
% 	 dx, dy, dg, ax, ay, ag, dx_alf, ag_alf, ag_msq;
% 	 p.cpg1, p.cpg1_msq, p.cpg1_alf, p.cpg2, p.cpg2_msq, p.cpg2_alf;
	 sa = sin(p.alfa);
	 ca = cos(p.alfa);

	p.xcp = 0.0;

	p.beta     = sqrt(1.0 - p.minf*p.minf);
	p.beta_msq = -0.5/p.beta;
	
	p.bfac	 = 0.5*p.minf*p.minf / (1.0 + p.beta);
	p.bfac_msq = 0.5/ (1.0 + p.beta)- p.bfac/ (1.0 + p.beta) * p.beta_msq;
	
	p.cl = 0.0;
	p.cm = 0.0;
	
	p.cdp = 0.0;
	
	p.cl_alf = 0.0;
	p.cl_msq = 0.0;
	
	i = 1;
	p.cginc = 1.0 - (p.gam(i)/p.qinf)*(p.gam(i)/p.qinf);
	p.cpg1	 = p.cginc/(p.beta + p.bfac*p.cginc);
	p.cpg1_msq = -p.cpg1/(p.beta + p.bfac*p.cginc)*(p.beta_msq + p.bfac_msq*p.cginc);
	
	p.cpi_gam = -2.0*p.gam(i)/p.qinf/p.qinf;
	p.cpc_cpi = (1.0 - p.bfac*p.cpg1)/ (p.beta + p.bfac*p.cginc);
	p.cpg1_alf = p.cpc_cpi*p.cpi_gam*p.gam_a(i);
	
	for (i=1:p.n)
		ip = i+1;
		if(i==p.n) 
            ip = 1;
        end
		
		p.cginc      = 1.0 - (p.gam(ip)/p.qinf)*(p.gam(ip)/p.qinf);
		p.cpg2	   = p.cginc/(p.beta + p.bfac*p.cginc);
		p.cpg2_msq   = -p.cpg2/(p.beta + p.bfac*p.cginc)*(p.beta_msq + p.bfac_msq*p.cginc);
		
		p.cpi_gam    = -2.0*p.gam(ip)/p.qinf/p.qinf;
		p.cpc_cpi    = (1.0 - p.bfac*p.cpg2)/ (p.beta + p.bfac*p.cginc);
		p.cpg2_alf   = p.cpc_cpi*p.cpi_gam*p.gam_a(ip);
		
		dx = (p.x(ip) - p.x(i))*ca + (p.y(ip) - p.y(i))*sa;
		dy = (p.y(ip) - p.y(i))*ca - (p.x(ip) - p.x(i))*sa;
		dg = p.cpg2 - p.cpg1;
		
		ax = (0.5*(p.x(ip)+p.x(i))-xref)*ca + (0.5*(p.y(ip)+p.y(i))-yref)*sa;
		ay = (0.5*(p.y(ip)+p.y(i))-yref)*ca - (0.5*(p.x(ip)+p.x(i))-xref)*sa;
		ag = 0.5*(p.cpg2 + p.cpg1);
		
		dx_alf = -(p.x(ip) - p.x(i))*sa + (p.y(ip) - p.y(i))*ca;
		ag_alf = 0.5*(p.cpg2_alf + p.cpg1_alf);
		ag_msq = 0.5*(p.cpg2_msq + p.cpg1_msq);
		
		p.cl	   = p.cl 	+ dx* ag;
		p.cdp    = p.cdp	- dy* ag;
		p.cm	   = p.cm 	- dx*(ag*ax + dg*dx/12.0)...
						- dy*(ag*ay + dg*dy/12.0);

		p.xcp = p.xcp+ dx*ag*(p.x(ip)+p.x(i))/2.0;
			
		p.cl_alf = p.cl_alf + dx*ag_alf + ag*dx_alf;
		p.cl_msq = p.cl_msq + dx*ag_msq;
		
		p.cpg1 = p.cpg2;
		p.cpg1_alf = p.cpg2_alf;
		p.cpg1_msq = p.cpg2_msq;
    end

	if(abs(p.cl)>0.0) 
        p.xcp= p.xcp/ p.cl;
    else
        p.xcp = 0.0;
    end
% 	return true;
end