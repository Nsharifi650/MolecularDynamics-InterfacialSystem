clear all
close all
%% shifting the values
T_shift = 2 ;% atom type shift - ALSO MASS
B_shift = 3 ;% bond 
A_shift = 1 ;% angle
D_shift = 0 ;%Dihedral
I_shift = 0 ;% Improper

%% import data files
Position=load('Moltemplate\position.dat');
angle=load('Moltemplate\angles.dat');
bond = load('Moltemplate\bonds.dat');
dihedrals = load('Moltemplate\dihedral.dat');
mass = load('Moltemplate\mass.dat');
impropers = load('Moltemplate\improper.dat');
charges = load('Moltemplate\Charge.dat');

%% import coefficient files
Coe_p=load('Coeffs/PairCoeff.dat');
Coe_a=load('Coeffs/angleCoeff.dat');
Coe_b = load('Coeffs/bondCoeff.dat');
Coe_d = load('Coeffs/DihedralCoeff.dat');
Coe_i = load('Coeffs/improperCoeff.dat');

%% position and paircoefficient

% positon
AtomTypeList = unique(Position(:,3),'stable');
for i= 1:length(Position(:,1))
    AtomType = Position(i,3);
    Cindex = find(charges==AtomType);
    Tindex = find(AtomTypeList==AtomType);
    Position(i,4) = charges(Cindex,2);
    Position(i,3) = Tindex;  
end
Position(:,3) = Position(:,3)+T_shift;

% Masses
for i= 1:length(AtomTypeList)
    m_index = find(mass(:,1)==AtomTypeList(i));
    MassOut(i,1)= i;
    MassOut(i,2)= mass(m_index,2);
end
MassOut(:,1) = MassOut(:,1) + T_shift;

% pair coefficient
for i= 1:length(AtomTypeList)
    P_C_index = find(Coe_p(:,1)==AtomTypeList(i));
    Pair_coeff(i,1:2)= i;
    Pair_coeff(i,3:4)= Coe_p(P_C_index,3:4);
end
Pair_coeff(:,1:2) = Pair_coeff(:,1:2)+ T_shift;

%% angles
% angle
AngleTypeList = unique(angle(:,2),'stable');
for i= 1:length(angle(:,1))
    AngleType = angle(i,2);
    Tindex = find(AngleTypeList==AngleType);
    angle(i,2) = Tindex;
end
Nangletypes = length(AngleTypeList);
angle(:,2) = angle(:,2) + A_shift;

% angle coefficient
for i= 1:length(AngleTypeList)
    A_C_index = find(Coe_a(:,1)==AngleTypeList(i));
    angle_coeff(i,1)= i;
    angle_coeff(i,2:3)= Coe_a(A_C_index,2:3);
end
angle_coeff(:,1)= angle_coeff(:,1)+A_shift;
%% bonds
% bond 
BondTypeList = unique(bond(:,2),'stable');
for i= 1:length(bond(:,1))
    bondType = bond(i,2);
    Tindex = find(BondTypeList==bondType);
    bond(i,2) = Tindex;
end
Nbondtypes = length(BondTypeList);
bond(:,2) = bond(:,2) + B_shift;

% bond coefficient
for i= 1:length(BondTypeList)
    B_C_index = find(Coe_b(:,1)==BondTypeList(i));
    bond_coeff(i,1)= i;
    bond_coeff(i,2:3)= Coe_b(B_C_index,2:3);
end
bond_coeff(:,1) = bond_coeff(:,1) + B_shift;

%% dihedrals
% dihedrals
DihedralTypeList = unique(dihedrals(:,2),'stable');
for i= 1:length(dihedrals(:,1))
    dihedralType = dihedrals(i,2);
    Tindex = find(DihedralTypeList==dihedralType);
    dihedrals(i,2) = Tindex;
end
Ndihedralstypes = length(DihedralTypeList);
dihedrals(:,2) = dihedrals(:,2) + D_shift;

%dihedral coefficient
for i= 1:length(DihedralTypeList)
    D_C_index = find(Coe_d(:,1)==DihedralTypeList(i));
    dihedral_coeff(i,1)= i;
    dihedral_coeff(i,2:5)= Coe_d(D_C_index,2:5);
end

dihedral_coeff(:,1) = dihedral_coeff(:,1) + D_shift;


%% impropers
% improper
ImproperTypeList = unique(impropers(:,2),'stable');

for i= 1:length(impropers(:,1))
    improperType = impropers(i,2);
    Tindex = find(ImproperTypeList==improperType);
    impropers(i,2) = Tindex;
end

% improper coefficient
for i= 1:length(ImproperTypeList)
    I_C_index = find(Coe_i(:,1)==ImproperTypeList(i));
    improper_coeff(i,1)= i;
    improper_coeff(i,2:4)= Coe_i(I_C_index,2:4);
end

%% write up molecule data files
Pfid = fopen('position.dat','wt');
for ii=1:length(Position(:,1))
	for jj=1:7
		fprintf(Pfid, num2str(Position(ii,jj)));
		fprintf(Pfid, '	');
	end
	fprintf(Pfid, '\n');
end
fclose(Pfid);

Bfid = fopen('Bond.dat','wt');
for ii=1:length(bond(:,1))
	for jj=1:4
		fprintf(Bfid, num2str(bond(ii,jj)));
		fprintf(Bfid, '	');
	end
	fprintf(Bfid, '\n');
end
fclose(Bfid);

Afid = fopen('Angle.dat','wt');
for ii=1:length(angle(:,1))
	for jj=1:5
		fprintf(Afid, num2str(angle(ii,jj)));
		fprintf(Afid, '	');
	end
	fprintf(Afid, '\n');
end
fclose(Afid);

Dfid = fopen('Dihedral.dat','wt');
for ii=1:length(dihedrals(:,1))
	for jj=1:6
		fprintf(Dfid, num2str(dihedrals(ii,jj)));
		fprintf(Dfid, '	');
	end
	fprintf(Dfid, '\n');
end
fclose(Dfid);

Mfid = fopen('Mass.dat','wt');
for ii=1:length(MassOut(:,1))
	for jj=1:2
		fprintf(Mfid, num2str(MassOut(ii,jj)));
		fprintf(Mfid, '	');
	end
	fprintf(Mfid, '\n');
end
fclose(Mfid);

ifid = fopen('Improper.dat','wt');
for ii=1:length(impropers(:,1))
	for jj=1:6
		fprintf(ifid, num2str(impropers(ii,jj)));
		fprintf(ifid, '	');
	end
	fprintf(ifid, '\n');
end
fclose(ifid);

%% writing the ceoff

Dfid = fopen('PARM.lammps','wt');
fprintf(Dfid, '#masses\n\n');
for ii=1:length(MassOut(:,1))
    fprintf(Dfid,"mass ");
	for jj=1:2
		fprintf(Dfid, num2str(MassOut(ii,jj)));
		fprintf(Dfid, '	');
	end
	fprintf(Dfid, '\n');
end
fprintf(Dfid, '\n');

fprintf(Dfid, '#Pair Coeffs\n\n');
for ii=1:length(Pair_coeff(:,1))
    fprintf(Dfid,"pair_coeff ");
	for jj=1:4
		fprintf(Dfid, num2str(Pair_coeff(ii,jj)));
		fprintf(Dfid, '	');
	end
	fprintf(Dfid, '\n');
end
fprintf(Dfid, '\n');

fprintf(Dfid, '#Bond Coeffs\n\n');
for ii=1:length(bond_coeff(:,1))
    fprintf(Dfid,"bond_coeff ");
	for jj=1:3
		fprintf(Dfid, num2str(bond_coeff(ii,jj)));
		fprintf(Dfid, '	');
	end
	fprintf(Dfid, '\n');
end
fprintf(Dfid, '\n');

fprintf(Dfid, '#Angle Coeffs\n\n');
for ii=1:length(angle_coeff(:,1))
    fprintf(Dfid,"angle_coeff ");
	for jj=1:3
		fprintf(Dfid, num2str(angle_coeff(ii,jj)));
		fprintf(Dfid, '	');
	end
	fprintf(Dfid, '\n');
end
fprintf(Dfid, '\n');

fprintf(Dfid, '#dihedral coeff\n\n');
for ii=1:length(dihedral_coeff(:,1))
    fprintf(Dfid,"dihedral_coeff ");
	for jj=1:5
		fprintf(Dfid, num2str(dihedral_coeff(ii,jj)));
		fprintf(Dfid, '	');
	end
	fprintf(Dfid, '\n');
end
fprintf(Dfid, '\n');

fprintf(Dfid, '#improper coeff\n\n');
for ii=1:length(improper_coeff(:,1))
    fprintf(Dfid,"improper_coeff ");
	for jj=1:4
		fprintf(Dfid, num2str(improper_coeff(ii,jj)));
		fprintf(Dfid, '	');
	end
	fprintf(Dfid, '\n');
end

disp('Done writing the parameter file')
fclose(Dfid);


