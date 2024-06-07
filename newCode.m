% Chain Defination
rhostar=0.0675;  % Density of Box
iseed=592984;    % No of seeds (8 digit or less)
nsets=1;         % No of sets of chain (blank line + 6 values for each sets)
swaptype=0;     % molecule tag rule: 0=by  mole, 1= from 1 end, 2= from 2nd end
nchain=40000;      % number of chains
nmonomer=1;    % monomer/chain
ntype=1;         % type of monomers (from output into LAMMPS file)
nbondtype=1;     % type of bonds (from output into LAMMPS file)
bondlength=1.0; % distance between monomers in reduced units
restrict=1.0;   % no distance from less than this from site i-1 to i+1 (reduced units)

natoms=0;
natoms=natoms+nchain*nmonomer;

%Set box size (sigma=1)
volume=natoms/rhostar;
xprd=volume^(1/3);
yprd=xprd;
zprd=xprd;
xboundlo = -xprd/2.0;
xboundhi = -xboundlo;
yboundlo = xboundlo;
yboundhi = xboundhi;
zboundlo = xboundlo;
zboundhi = xboundhi;

%generate random chains
% loop over set and chain in each set
n=0;
nmolecule=0;
    for ichain=1:nchain
        nmolecule= nmolecule+1;      

 % random starting point for the chain in the box
        x1=0;
        y1=0;
        z1=0;
        x2= xboundlo + rand(1)*xprd;
        y2= yboundlo + rand(1)*yprd;
        z2= zboundlo + rand(1)*zprd;
        
% store 1st monomer of chain
% 1st monomer is always in original box (image=0)
    [x2,y2,z2]=pbc(x2,y2,z2,xboundlo,yboundlo,zboundlo,xboundhi,yboundhi,zboundhi,xprd,yprd,zprd);
       n=n+1;
       x(n)=round(x2*10^4)/10^4;
       y(n)=round(y2*10^4)/10^4;
       z(n)=round(z2*10^4)/10^4;

          imagex(n) = 0;
          imagey(n) = 0;
          imagez(n) = 0;
          if (swaptype == 0)
            molecule(n) = nmolecule ;
          else
            molecule(n) = 1 ;
          end
   % generate rest of monomers in this chain 
  for imonomer=2:nmonomer

            x0 = x1;
            y0 = y1;
            z0 = z1;
            x1 = x2;
            y1 = y2;
            z1 = z2;
  

   % random point inside sphere of unit radius  
   while true 
   while true   
            xinner = 2.0*rand(1) - 1.0;
            yinner = 2.0*rand(1) - 1.0;
            zinner = 2.0*rand(1) - 1.0;  
            rsq = xinner * xinner + yinner * yinner + zinner * zinner;
   
   if(rsq<=1)
   break; 
   end
   end
% project point to surface of sphere of unit radius
            r = (rsq)^(1/2) ;
            xsurf = xinner/r;
            ysurf = yinner/r;
            zsurf = zinner/r;

% create new point by scaling unit offsets by bondlength (sigma = 1.0)
            x2 = x1 + round(xsurf*bondlength*10^4)/10^4;
            y2 = y1 + round(ysurf*bondlength*10^4)/10^4;
            z2 = z1 + round(zsurf*bondlength*10^4)/10^4;

%check that new point meets restriction requirement only for 3rd monomer and beyond 

            dx = x2 - x0;
            dy = y2 - y0;
            dz = z2 - z0;
             r = sqrt(dx*dx + dy*dy + dz*dz);

   if((imonomer<=2) || (r > restrict))
   break;
   end
   end
% store new points
%  if delta to previous bead is large, then increment/decrement image flag
[x2,y2,z2]=pbc(x2,y2,z2,xboundlo,yboundlo,zboundlo,xboundhi,yboundhi,zboundhi,xprd,yprd,zprd);
       n=n+1;
       x(n)=round(x2*10^4)/10^4;
       y(n)=round(y2*10^4)/10^4;
       z(n)=round(z2*10^4)/10^4;
     if (abs(x(n)-x(n-1)) < 2.0*bondlength)
              imagex(n) = imagex(n-1);
            elseif (x(n) - x(n-1) < 0.0)
              imagex(n) = imagex(n-1) + 1;
            elseif (x(n) - x(n-1) > 0.0)
              imagex(n) = imagex(n-1) - 1;
     end
    if (abs(y(n)-y(n-1)) < 2.0*bondlength)
              imagey(n) = imagey(n-1);
            elseif (y(n) - y(n-1) < 0.0)
              imagey(n) = imagey(n-1) + 1;
            elseif (y(n) - y(n-1) > 0.0)
              imagey(n) = imagey(n-1) - 1;
    end
             
            if (abs(z(n)-z(n-1)) < 2.0*bondlength)
              imagez(n) = imagez(n-1) ;
            elseif (z(n) - z(n-1) < 0.0)
              imagez(n) = imagez(n-1) + 1;
            elseif (z(n) - z(n-1) > 0.0)
              imagez(n) = imagez(n-1) - 1;
            end

if (swaptype == 0)
              molecule(n) = nmolecule;
            elseif (swaptype == 1)
              molecule(n) = imonomer;
            elseif (swaptype == 2)
              if (imonomer <= nmonomer/2)
                molecule(n) = imonomer;
              else
                molecule(n) = nmonomer+1-imonomer;
              end
end
end
end
% c compute quantities needed for LAMMPS file
nbonds = 0;
ntypes = 0;
nbondtypes = 0;
for i=1:1
nbonds = nbonds + nchain*(nmonomer-1);
        if (ntype > ntypes) 
            ntypes = ntype;
        if (nbondtype > nbondtypes)
            nbondtypes = nbondtype;
end
end
end
nangle=0;
ndihedrals=0;
nimpropers=0;
nangletypes=0;
ndihedralstypes=0;
nimpropertypes=0;
%Atoms matrix
for i=1:natoms
E(i,1)=i;
E(i,2)=molecule(i);
E(i,3)=1;
E(i,4)=x(i);
E(i,5)=y(i);
E(i,6)=z(i);
E(i,7)=imagex(i);
E(i,8)=imagey(i);
E(i,9)=imagez(i);
end

% write out LAMMPS file
FID = fopen('Fortran_LowMoleculedata.txt', 'w');
fprintf(FID, 'Model for Poly Ethylene\n');
fprintf(FID, '\n');
fprintf(FID, '%g\tatoms\n',natoms);
fprintf(FID, '%g\tbonds\n',nbonds);
fprintf(FID, '%g\tangles\n',nangle);
fprintf(FID, '%g\tdihedrals\n',ndihedrals);
fprintf(FID, '%g\timpropers\n',nimpropers);
fprintf(FID, '\n');
fprintf(FID, '%g\tatom types\n',ntypes);
fprintf(FID, '%g\tbond types\n',nbondtypes);
fprintf(FID, '%g\tangle types\n',nangletypes);
fprintf(FID, '%g\tdihedral types\n',ndihedralstypes);
fprintf(FID, '%g\timproper types\n',nimpropertypes);
fprintf(FID, '\n');
fprintf(FID, '%g\t%g\txlo\txhi\n',xboundlo,xboundhi);
fprintf(FID, '%g\t%g\tylo\tyhi\n',yboundlo,yboundhi);
fprintf(FID, '%g\t%g\tzlo\tzhi\n',zboundlo,zboundhi);
fprintf(FID, '\n');
fprintf(FID, 'Masses\n');
fprintf(FID, '\n');
fprintf(FID, '1\t1.0\n');
fprintf(FID, '\n');
fprintf(FID, 'Atoms\n');
fprintf(FID, '\n');
if FID == -1, error('Cannot create file.'); end
fprintf(FID, '%g %g %g %g %g %g %g %g %g\n', E');
fprintf(FID, '\n');
fprintf(FID, 'Bonds\n');
fprintf(FID, '\n');
%Bond Formation
n=0;
m=0;
for ichain=1:nchain
for imonomer=1:nmonomer
n=n+1;
if(imonomer~=nmonomer)
    m=m+1;
fprintf(FID, '%g %g %g %g\n',m,nbondtype,n,n+1);
end
end
end
fclose(FID);