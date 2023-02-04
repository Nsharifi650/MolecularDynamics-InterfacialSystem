clear all
close all

%% setting up the walls 
% rotating crystal
X_rot = 0;
Y_rot = 0;
Z_rot = 0;
Surface=importdata('MS_files_surface/Diamond.dat');
Surface = Surface.textdata;
Sur_type = Surface(:,8);
Sur_pos = str2double(Surface(:,2:4));

syms t
Rx = [1 0 0; 0 cos(t) -sin(t); 0 sin(t) cos(t)];
Ry = [cos(t) 0 sin(t); 0 1 0; -sin(t) 0 cos(t)];
Rz = [cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1];
for jj=1:length(Surface(:,2))
    Sur_rot1 = Rx*[Sur_pos(jj,1);Sur_pos(jj,2);Sur_pos(jj,3)];
    Rotation1 = subs(Sur_rot1, t, X_rot);
    Rotation2 = subs(Ry*Rotation1, t, Y_rot);
    Rotation3 = subs(Rz*Rotation2, t, Z_rot);
    Sur_pos(jj,:) = Rotation3;
end

% recentre
x_shift = max(Sur_pos(:,1))-min(Sur_pos(:,1));
y_shift = max(Sur_pos(:,2))-min(Sur_pos(:,2));
z_shift = max(Sur_pos(:,3))-min(Sur_pos(:,3));
Sur_pos(:,1) = Sur_pos(:,1) - min(Sur_pos(:,1))-x_shift/2;
Sur_pos(:,2) = Sur_pos(:,2) - min(Sur_pos(:,2))-y_shift/2;
Sur_pos(:,3) = Sur_pos(:,3) - min(Sur_pos(:,3)); % moving bottom layer to zeros


%% system size
dFeO=2.085;
widthwall=6*dFeO;
d_gap= 150;
walldistance=d_gap;

% graphene
c=20;

% water
ee=6;

%% top wall 
% position of lowest carbon atom
Bottom_atom_position = min(Sur_pos(:,3));
id_atom = 0;
id_mol = 1;
xx= 1;
for ii=1:length(Surface(:,2))
    id_atom=id_atom+1;
    if strcmp(Sur_type(ii),'C')
        if Sur_pos(ii,3) <= Bottom_atom_position+0.5
            id_type = 20;
            charge = 0.2650;
            Bottom_layer_index(xx) = ii;
            xx = xx+1;
        else
        id_type = 21;
        charge = 0;
        end
    end
    A(ii,:)=[id_atom id_mol id_type charge Sur_pos(ii,1) Sur_pos(ii,2) Sur_pos(ii,3)+(d_gap)/2];
end

A_topwall = A; % top wall before adding OH layer
% adding the OH layer to the top wall:
O_charge = -0.683;
H_charge = 0.418;
C_num = length(A(:,1));
xx = 1;
for ii=1:length(Surface(:,2))
    if A(ii,3) == 20
        C_layer(xx,:)= A(ii,:);
        O_layer(xx,:) = [C_num+xx id_mol 22 O_charge A(ii,5) A(ii,6) A(ii,7)-1];
        H_layer(xx,:) = [C_num+length(Bottom_layer_index)+xx id_mol 23 H_charge A(ii,5) A(ii,6) A(ii,7)-2];
        xx=xx+1;
    end
end

A = [A ; O_layer; H_layer];

% wall 2
id_atom = length(A(:,1));
id_mol = id_mol +1;

At=A;
At(:,7)=-At(:,7);
At(:,2)=id_mol;
At(:,1)=At(:,1)+id_atom;
A=[A; At];
% A(:,3)=A(:,3)+2; %% shifting atom types !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%% creating bonds for the Diamond
%n is the number of atoms in the solid
n= length(A_topwall(:,1));
R_min = 4;
xx = 1;
for i=1:n
    p1 = A(i,5:7); % position of atom of interest
    p_rest = A(1:n,5:7) - p1; % position of the rest of atoms relative to p1
    radius = sqrt(p_rest(:,1).^2+p_rest(:,2).^2+p_rest(:,3).^2);
    for j = 1:length(radius)
        if radius(j) < R_min
            if i == j
            else
                S_B(xx,:) = [xx xx i j];
                P_bond(xx,:) = [xx 2000 radius(j)];
                xx= xx+1;
            end
        end
    end
end
% including the second wall
S_B_2ndwall = S_B;
mm = length(S_B(:,1));
S_B_2ndwall(:,1:2) = S_B_2ndwall(:,1:2)+mm; % shifting bond types
S_B_2ndwall(:,3:4) = S_B_2ndwall(:,3:4)+n + 2*length(O_layer); % shifting atom numbers


P_bond_2ndwall = P_bond;
P_bond_2ndwall(:,1) = P_bond_2ndwall(:,1)+mm;

S_B = [S_B; S_B_2ndwall];
P_bond = [P_bond; P_bond_2ndwall];

%% creating COH bonds at the surface

mm = length(O_layer);
for i =1:mm
    S_CO_bonds(i,:) =  [i 1 C_layer(i,1) O_layer(i,1)]; % bond num, type, atom i, atom j
    S_OH_bonds(i,:) =  [i 1 O_layer(i,1) H_layer(i,1)];
end

% OH layer of TOP wall: shifting bond numbers and bond_type before adding to the bond list
Bond_Numbershift = length(S_B(:,1));
Bond_typeshift = max(S_B(:,2));
S_CO_bonds(:,1)= S_CO_bonds(:,1)+Bond_Numbershift;
S_CO_bonds(:,2)= S_CO_bonds(:,2)+Bond_typeshift;

S_OH_bonds(:,1)= S_OH_bonds(:,1)+Bond_Numbershift+length(S_CO_bonds(:,1));
S_OH_bonds(:,2)= S_OH_bonds(:,2)+Bond_typeshift+1;

S_B=[S_B; S_CO_bonds; S_OH_bonds];

%Parameters file

CO_para = 320;
CO_bond_length = 1.41;
OH_para = 553;
OH_bond_length = 0.945;

P_CO_bond = [Bond_typeshift+1 CO_para CO_bond_length];
P_OH_bond = [Bond_typeshift+2 OH_para OH_bond_length];

P_bond = [P_bond; P_CO_bond; P_OH_bond];


% OH layer of Bottom wall: shifting bond numbers and bond_type before adding to the bond list
aa = length(A_topwall(:,1)) + 2*length(O_layer);

mm = length(O_layer);

for i =1:mm
    S_CO_bonds(i,:) =  [i 1 C_layer(i,1)+aa O_layer(i,1)+aa]; % bond num, type, atom i, atom j
    S_OH_bonds(i,:) =  [i 1 O_layer(i,1)+aa H_layer(i,1)+aa];
end

Bond_Numbershift = length(S_B(:,1));
Bond_typeshift = max(S_B(:,2));
S_CO_bonds(:,1)= S_CO_bonds(:,1)+Bond_Numbershift;
S_CO_bonds(:,2)= S_CO_bonds(:,2)+Bond_typeshift;

S_OH_bonds(:,1)= S_OH_bonds(:,1)+Bond_Numbershift+length(S_CO_bonds(:,1));
S_OH_bonds(:,2)= S_OH_bonds(:,2)+Bond_typeshift+1;

S_B=[S_B; S_CO_bonds; S_OH_bonds];

%Parameters file

CO_para = 320;
CO_bond_length = 1.41;
OH_para = 553;
OH_bond_length = 0.945;

P_CO_bond = [Bond_typeshift+1 CO_para CO_bond_length];
P_OH_bond = [Bond_typeshift+2 OH_para OH_bond_length];

P_bond = [P_bond; P_CO_bond; P_OH_bond];

txlo=min(A(:,5)); txhi=-txlo;
tylo=min(A(:,6)); tyhi=-tylo;
tzlo=min(A(:,7)); tzhi=-tzlo;


%% solutions

x=txlo+dFeO/2;
y=tylo+dFeO/2;
z=-widthwall/2;

% load  additive
cptatom=max(A(:,1));
cptmol= id_mol;

cptbond=0; 
cptangle=0; 
cptdihedrals=0; 
cptimpropers=0; 

PnmpA=load('AOT/position.dat');
PnmpA(:,5)=PnmpA(:,5)-mean(PnmpA(:,5));
PnmpA(:,6)=PnmpA(:,6)-mean(PnmpA(:,6));
PnmpA(:,7)=PnmpA(:,7)-mean(PnmpA(:,7));
BnmpA=load('AOT/Bond.dat');
AnmpA=load('AOT/Angle.dat');
DnmpA=load('AOT/Dihedral.dat');
InmpA=load('AOT/Improper.dat');
%%%% rotating the dodecanoic acid molecule

% X_rot = pi/2;
% Y_rot = -pi/2+(-pi/2);
% Z_rot = pi;


X_rot = pi/2;
Y_rot = -pi/2+(-pi/2);
Z_rot = pi;

syms t
Rx = [1 0 0; 0 cos(t) -sin(t); 0 sin(t) cos(t)];
Ry = [cos(t) 0 sin(t); 0 1 0; -sin(t) 0 cos(t)];
Rz = [cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1];
for jj=1:length(PnmpA(:,5))
    Sur_rot1 = Rx*[PnmpA(jj,5);PnmpA(jj,6);PnmpA(jj,7)];
    Rotation1 = subs(Sur_rot1, t, X_rot);
    Rotation2 = subs(Ry*Rotation1, t, Y_rot);
    Rotation3 = subs(Rz*Rotation2, t, Z_rot);
    PnmpA(jj,5:7) = Rotation3;
end

%%%%%%
z = -60;
y = tylo+3;
x = txlo+5;
N_bilayers = 2;
for mm=1:2
    for ii=1:9
        for jj=1:3
            cptatom0=cptatom;
            cptmol=cptmol+1;
            for ii=1:length(BnmpA(:,1))
	            cptbond=cptbond+1;
	            B(cptbond,:)=[cptbond BnmpA(ii,2) BnmpA(ii,3)+cptatom BnmpA(ii,4)+cptatom];
            end
            for ii=1:length(AnmpA(:,1))
	            cptangle=cptangle+1;
	            Ag(cptangle,:)=[cptangle AnmpA(ii,2) AnmpA(ii,3)+cptatom AnmpA(ii,4)+cptatom AnmpA(ii,5)+cptatom];
            end
            for ii=1:length(DnmpA(:,1))
	            cptdihedrals=cptdihedrals+1;
	            D(cptdihedrals,:)=[cptdihedrals DnmpA(ii,2) DnmpA(ii,3)+cptatom DnmpA(ii,4)+cptatom DnmpA(ii,5)+cptatom DnmpA(ii,6)+cptatom];
            end
            
            for ii=1:length(InmpA(:,1))
                cptimpropers=cptimpropers+1;
                Ig(cptimpropers,:)=[cptimpropers InmpA(ii,2) InmpA(ii,3)+cptatom InmpA(ii,4)+cptatom InmpA(ii,5)+cptatom InmpA(ii,6)+cptatom];
            end
        
            for ii=1:length(PnmpA(:,1))
	            cptatom=cptatom+1;
	            A(cptatom,:)=[PnmpA(ii,1)+cptatom0 cptmol PnmpA(ii,3) PnmpA(ii,4) PnmpA(ii,5)+x PnmpA(ii,6)+y PnmpA(ii,7)+z];
            end
            x = x +16;
        end
        y = y+5;
        x = txlo+5;
    end
    z = z + 105;
    y = tylo+3;
    x = txlo+5;
end


%%%%%%% second part of the bilayer
% X_rot = -pi/2;
% Y_rot = -pi/2;
% Z_rot = pi;

X_rot = pi;
Y_rot = 0;
Z_rot = 0;
syms t
Rx = [1 0 0; 0 cos(t) -sin(t); 0 sin(t) cos(t)];
Ry = [cos(t) 0 sin(t); 0 1 0; -sin(t) 0 cos(t)];
Rz = [cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1];
for jj=1:length(PnmpA(:,5))
    Sur_rot1 = Rx*[PnmpA(jj,5);PnmpA(jj,6);PnmpA(jj,7)];
    Rotation1 = subs(Sur_rot1, t, X_rot);
    Rotation2 = subs(Ry*Rotation1, t, Y_rot);
    Rotation3 = subs(Rz*Rotation2, t, Z_rot);
    PnmpA(jj,5:7) = Rotation3;
end
% 
% %%%%%%%%%
% 
z = 60;
y = tylo+3;
x = txlo+5;
for mm=1:2
    for ii=1:9
        for jj=1:3
            cptatom0=cptatom;
            cptmol=cptmol+1;
            for ii=1:length(BnmpA(:,1))
	            cptbond=cptbond+1;
	            B(cptbond,:)=[cptbond BnmpA(ii,2) BnmpA(ii,3)+cptatom BnmpA(ii,4)+cptatom];
            end
            for ii=1:length(AnmpA(:,1))
	            cptangle=cptangle+1;
	            Ag(cptangle,:)=[cptangle AnmpA(ii,2) AnmpA(ii,3)+cptatom AnmpA(ii,4)+cptatom AnmpA(ii,5)+cptatom];
            end
            for ii=1:length(DnmpA(:,1))
	            cptdihedrals=cptdihedrals+1;
	            D(cptdihedrals,:)=[cptdihedrals DnmpA(ii,2) DnmpA(ii,3)+cptatom DnmpA(ii,4)+cptatom DnmpA(ii,5)+cptatom DnmpA(ii,6)+cptatom];
            end
            
            for ii=1:length(InmpA(:,1))
                cptimpropers=cptimpropers+1;
                Ig(cptimpropers,:)=[cptimpropers InmpA(ii,2) InmpA(ii,3)+cptatom InmpA(ii,4)+cptatom InmpA(ii,5)+cptatom InmpA(ii,6)+cptatom];
            end
        
            for ii=1:length(PnmpA(:,1))
	            cptatom=cptatom+1;
	            A(cptatom,:)=[PnmpA(ii,1)+cptatom0 cptmol PnmpA(ii,3) PnmpA(ii,4) PnmpA(ii,5)+x PnmpA(ii,6)+y PnmpA(ii,7)+z];
            end
            x = x + 16;
        end
        y = y+5;
        x = txlo+5;
    end
    z = z -105;
    y = tylo+5;
    x = txlo+5;
end

% load Solvent
lengthAbeforewater=length(A(:,1));
clear PnmpA;
clear BnmpA;
clear AnmpA;
PnmpA=load('water/position.dat');
PnmpA(:,5)=PnmpA(:,5)-mean(PnmpA(:,5));
PnmpA(:,6)=PnmpA(:,6)-mean(PnmpA(:,6));
PnmpA(:,7)=PnmpA(:,7)-mean(PnmpA(:,7));
BnmpA=load('water/Bond.dat');
AnmpA=load('water/Angle.dat');
%DnmpA=load('Cyclohexane/Dihedral.dat');
%InmpA=load('./Hexadecane/Impropers.dat');
dx=0.4*ee; dy=0.4*ee; dz=0.4*ee;
x=txlo+dx;
y=tylo+dy;
z=-45;
while z<45 
	while y < tyhi-dy
		while x < txhi-dx
			xeff=x+(rand-0.5)/10;
			yeff=y+(rand-0.5)/10;
			zeff=z+(rand-0.5)/10;
		
			for ii=1:length(PnmpA(:,1))
				Glyc(ii,:)=[cptatom cptmol PnmpA(ii,3) PnmpA(ii,4) PnmpA(ii,5)+xeff PnmpA(ii,6)+yeff PnmpA(ii,7)+zeff];
			end
			mind=1e5;
			for ii=1:length(Glyc(:,1))
				for jj=1:lengthAbeforewater
					d=sqrt((Glyc(ii,5)-A(jj,5))^2+(Glyc(ii,6)-A(jj,6))^2+(Glyc(ii,7)-A(jj,7))^2);
					if d<mind
						mind=d;
					end
				end
			end
			if mind>3.4
				cptmol=cptmol+1;
				for ii=1:length(BnmpA(:,1))
					cptbond=cptbond+1;
					B(cptbond,:)=[cptbond BnmpA(ii,2) BnmpA(ii,3)+cptatom BnmpA(ii,4)+cptatom];
				end

				for ii=1:length(AnmpA(:,1))
					cptangle=cptangle+1;
					Ag(cptangle,:)=[cptangle AnmpA(ii,2) AnmpA(ii,3)+cptatom AnmpA(ii,4)+cptatom AnmpA(ii,5)+cptatom];
                end
%                 for ii=1:length(DnmpA(:,1))
%                     cptdihedrals=cptdihedrals+1;
%                     D(cptdihedrals,:)=[cptdihedrals DnmpA(ii,2) DnmpA(ii,3)+cptatom DnmpA(ii,4)+cptatom DnmpA(ii,5)+cptatom DnmpA(ii,6)+cptatom];
%                 end
%                 for ii=1:length(InmpA(:,1))
%                     cptimpropers=cptimpropers+1;
%                     Ig(cptimpropers,:)=[cptimpropers InmpA(ii,2) InmpA(ii,3)+cptatom InmpA(ii,4)+cptatom InmpA(ii,5)+cptatom InmpA(ii,6)+cptatom];
%                 end
                for ii=1:length(PnmpA(:,1))
					cptatom=cptatom+1;
					A(cptatom,:)=[cptatom cptmol PnmpA(ii,3) PnmpA(ii,4) PnmpA(ii,5)+xeff PnmpA(ii,6)+yeff PnmpA(ii,7)+zeff];
                end
			end
			x=x+dx;
		end
		x=txlo+dx/2;
		y=y+dy;
	end
	y=tylo+dy/2;
	z=z+dz;	
end
%%
% % % shifting the bond numbers of the solid
Bond_Numbershift = length(B(:,1));
Bond_typeshift = max(B(:,2));
S_B(:,1)= S_B(:,1)+Bond_Numbershift;
S_B(:,2)= S_B(:,2)+Bond_typeshift;
B=[B; S_B];

%shifting the solid param files
P_bond(:,1) = P_bond(:,1)+Bond_typeshift;


Natomtypes=max(A(:,3));
Nbondtypes=max(B(:,2)); %6;
Nangletypes=max(Ag(:,2)); %9;
Ndihedraltypes=max(D(:,2)); %10;  
Nimpropertypes=1;
cptbond = length(B(:,1));


% writing bond parm file for solid
Dfid = fopen('solidBonds.lammps','wt');

fprintf(Dfid, '#Bond Coeffs\n\n');
for ii=1:length(P_bond(:,1))
    fprintf(Dfid,"bond_coeff ");
	for jj=1:3
		fprintf(Dfid, num2str(P_bond(ii,jj)));
		fprintf(Dfid, '	');
	end
	fprintf(Dfid, '\n');
end
fprintf(Dfid, '\n');
fclose(Dfid);


fid = fopen('data.lammps','wt');
fprintf(fid, '# System\n\n');
fprintf(fid, num2str(cptatom));
fprintf(fid, ' atoms\n');
fprintf(fid, num2str(cptbond));
fprintf(fid, ' bonds\n');
fprintf(fid, num2str(cptangle));
fprintf(fid, ' angles\n');
fprintf(fid, num2str(cptdihedrals));
fprintf(fid, ' dihedrals\n');
fprintf(fid, num2str(cptimpropers));
fprintf(fid, ' impropers\n\n');
fprintf(fid, num2str(Natomtypes));
fprintf(fid, ' atom types\n');
fprintf(fid, num2str(Nbondtypes));
fprintf(fid, ' bond types\n');
fprintf(fid, num2str(Nangletypes));
fprintf(fid, ' angle types\n');
fprintf(fid, num2str(Ndihedraltypes));
fprintf(fid, ' dihedral types\n');
fprintf(fid, num2str(Nimpropertypes));
fprintf(fid, ' improper types\n');
fprintf(fid, 'extra bond per atom ');
fprintf(fid, num2str(2));
fprintf(fid, '\n');
fprintf(fid, 'extra angle per atom ');
fprintf(fid, num2str(1));
fprintf(fid, '\n');
fprintf(fid, 'extra special per atom ');
fprintf(fid, num2str(2));
fprintf(fid, '\n');
fprintf(fid, num2str([txlo txhi]));
fprintf(fid, ' xlo xhi\n');
fprintf(fid, num2str([tylo tyhi]));
fprintf(fid, ' ylo yhi\n');
fprintf(fid, num2str([tzlo tzhi]));
fprintf(fid, ' zlo zhi\n\n');
fprintf(fid, 'Atoms\n\n');
for ii=1:length(A(:,1))
	for jj=1:7
		fprintf(fid, num2str(A(ii,jj)));
		fprintf(fid, '	');
	end
	fprintf(fid, '\n');
end
fprintf(fid, '\n');
fprintf(fid, 'Bonds\n\n');
for ii=1:length(B(:,1))
	for jj=1:4
		fprintf(fid, num2str(B(ii,jj)));
		fprintf(fid, '	');
	end
	fprintf(fid, '\n');
end
fprintf(fid, '\n');
fprintf(fid, 'Angles\n\n');
for ii=1:length(Ag(:,1))
	for jj=1:5
		fprintf(fid, num2str(Ag(ii,jj)));
		fprintf(fid, '	');
	end
	fprintf(fid, '\n');
end
fprintf(fid, '\n');
fprintf(fid, 'Dihedrals\n\n');
for ii=1:length(D(:,1))
	for jj=1:6
		fprintf(fid, num2str(D(ii,jj)));
		fprintf(fid, '	');
	end
	fprintf(fid, '\n');
end
fprintf(fid, '\n');
fprintf(fid, 'Impropers\n\n');
for ii=1:length(Ig(:,1))
	for jj=1:6
		fprintf(fid, num2str(Ig(ii,jj)));
		fprintf(fid, '	');
	end
	fprintf(fid, '\n');
end
disp('Done writing the data file')
fclose(fid);