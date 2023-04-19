clear variables
close all

%inputs

filen = 'disp3_PHBST.nc'; %Anaddb output file
plotlengths = [24, 23, 21, 28, 16, 33, 22, 28, 16, 33, 23, 24]; %lenghts of segments in q space
plotlabels = {'Q_3', 'Q_2',  'Q_1', '\Gamma', 'M', 'K', '\Gamma', 'A', 'L', 'H', 'A|L', 'M|H', 'K'};
breakval = [69, 146, 245, 268]; %discontinuities in the q path
gammaval = [69, 146]; %Gamma points
plotflag = 1; %Set to 1 to produce a plot

%load stuff

eigend_abi = ncread(filen, 'phdispl_cart'); %Eigendisplacements
qpoints = ncread(filen, 'qpoints'); %List of qpoints
mass_amu = ncread(filen, 'atomic_mass_units'); %atomic masses in amu
species = ncread(filen, 'atom_species'); %atomic species labels
phf_eV = ncread(filen, 'phfreqs'); %Phonon energies in eV
bec = ncread(filen, 'becs_cart'); %Born effective charge tensors
qpts = ncread(filen, 'qpoints'); %Q points

pm = 1822.89; % 1 amu in electron mass units
b2A = 1.88973; % 1 Angstrom in Bohr radii
a2bohr = 0.529177; %1 Angstrom in Bohr radii
T2eV = 0.00413566972; % 1 THz in eV
hbar = 1; % Hartree units
hbar_SI = 1.054571817e-34; %J/Hz
e2C = 1.60218e-19; %1 e in C
emu2k = 9.10938e-31; %1 emu in kg
nmag_SI = 5.050783699e-27; %nuclear magneton in J/T
nmag = nmag_SI/1.85480201566e-23; %1 nuclear magneton in 
bmag_SI = 9.2740100783e-24; %Bohr magneton in J/T
bmag = bmag_SI/1.85480201566e-23; %Bohr magneton in Hartree units
eV2H = 0.037; %1 eV in Ha
Heb2V = 5.14220674763e11; %1 Ha/e/bohr in V/m
H2s = 2.4188843265857e-17; % convert atomic time units to s
eV2Hz = 241799050402293; % convert eV to Hz
permeability = 4*pi()*(1/137)^2;

natom = size(species,1); %number of atoms
nqpt = size(eigend_abi,4); %number of q points
nphonon = size(eigend_abi,2); %number of phonon bands
nspecies = size(mass_amu,1); %number of types of atoms

rprim = ncread(filen, 'primitive_vectors'); 
xred = transpose(ncread(filen,'reduced_atom_positions'));

vol = dot(rprim(:,1),cross(rprim(:,2),rprim(:,3)));

mass = mass_amu.*pm; %convert masses to emu

%Make xcart

for atom = 1:natom
    xcart(atom,:) = xred(atom,1)*rprim(1,:)+xred(atom,2)*rprim(2,:)+xred(atom,3)*rprim(3,:);
end

%Fix LO-TO splitting at gamma

if plotflag == 1

ngamma = size(gammaval,2);
nbreakval = size(breakval,2);

if ngamma >= 1
    
    phf_na = transpose(ncread(filen, 'non_analytical_phonon_modes'));
    eigend_abi_na = ncread(filen, 'non_analytical_phdispl_cart');
    isdone = zeros(size(gammaval));
    
    for gamma = 1:ngamma
        if gammaval(gamma) == 1
            phf_eV(:,gammaval(gamma)) = phf_na(gamma,:);
            eigend_abi(:,:,:,gammaval(gamma)) = eigend_abi_na(:,:,:,gamma);
            isdone(gamma) = 1;
        end
    end
    for gamma = 1:ngamma
        if gammaval(gamma) == nqpt
            phf_eV(:,gammaval(gamma)) = phf_na(ngamma,:);
            eigend_abi(:,:,:,gammaval(gamma)) = eigend_abi_na(:,:,:,ngamma);
            isdone(gamma) = 1;
        end
    end
    
    j = 0;
    
    for gamma = 1:ngamma
        if isdone(gamma) == 0
            phf_unbroken = phf_eV;
            phf_eV(:,gammaval(gamma)) = phf_na(gamma+j,:);
            phf_eV(:,gammaval(gamma)+1) = phf_na(gamma+j+1,:);
            
            for i = (gammaval(gamma)+2):size(phf_eV,2)+1
                phf_eV(:,i)= phf_unbroken(:,i-1);
            end
            
            eigend_unbroken = eigend_abi;
            eigend_abi(:,:,:,gammaval(gamma)) = eigend_abi_na(:,:,:,gamma+j);
            eigend_abi(:,:,:,gammaval(gamma)+1) = eigend_abi_na(:,:,:,gamma+j+1);
            
            for i = (gammaval(gamma)+2):size(eigend_abi,4)+1
                eigend_abi(:,:,:,i) = eigend_unbroken(:,:,:,i-1);
            end
            
            length = 0;
            
            for i = 1:size(plotlengths,2)
                length = length + plotlengths(1,i);
                if length == gammaval(gamma) - 1
                    plotlengths(1,i+1) = plotlengths(1,i+1) + 1;
                end
            end
            
            if nbreakval >= 1
                for i = 1:size(breakval,2)
                    if breakval(i) > gammaval(gamma)
                        breakval(i) = breakval(i) + 1;
                    end
                end
            end
            
            if ngamma >= 2
                for i = 1:size(gammaval,2)
                    if gammaval(i) > gammaval(gamma)
                        gammaval(i) = gammaval(i) + 1;
                    end
                end
            end
            
            nqpt = nqpt + 1;
            j = j + 1;
            isdone(gamma) = 1;
        end
    end
end

end
            
phf = transpose(phf_eV*(1000)); %convert to meV

%make R L and Z vectors

for atom = 1:natom
	R = zeros(1,natom*3);
	R(1,(atom-1)*3+1)= 1;
	R(1,(atom-1)*3+2) = 1i;
	Rz{atom} = 1*(2^-0.5).*R;
	Lz{atom} = 1*(2^-0.5).*conj(R);
	Z = zeros(1,natom*3);
	Z(1,(atom-1)*3+3) = 1;
	Zz{atom} = Z;
end

for atom = 1:natom
	R = zeros(1,natom*3);
	R(1,(atom-1)*3+3)= 1;
	R(1,(atom-1)*3+1) = 1i;
	Ry{atom} = 1*(2^-0.5).*R;
	Ly{atom} = 1*(2^-0.5).*conj(R);
	Z = zeros(1,natom*3);
	 Z(1,(atom-1)*3+2) = 1;
	Zy{atom} = Z;
end

for atom = 1:natom
	R = zeros(1,natom*3);
	R(1,(atom-1)*3+2)= 1;
	R(1,(atom-1)*3+3) = 1i;
	Rx{atom} = 1*(2^-0.5).*R;
	Lx{atom} = 1*(2^-0.5).*conj(R);
	Z = zeros(1,natom*3);
	Z(1,(atom-1)*3+1) = 1;
	Zx{atom} = Z;
end

%load eigenvectors (Angstrom*√emu) and calculate phonon circular polarization

for phonon = 1:nphonon
    for qpt = 1:nqpt
		for atom = 1:natom
			amf(atom) = sqrt((mass(species(atom))));
			for dir = 1:3
                eigendr(atom,dir,phonon,qpt) = eigend_abi(1,(atom-1)*3+dir,phonon,qpt)/a2bohr; %real part (bohr)
                eigendi(atom,dir,phonon,qpt) = eigend_abi(2,(atom-1)*3+dir,phonon,qpt)/a2bohr; %imaginary part
                eigend(atom,dir,phonon,qpt) = complex(eigendr(atom,dir,phonon,qpt),eigendi(atom,dir,phonon,qpt));
                eivecr(atom,dir,phonon,qpt) = eigendr(atom,dir,phonon,qpt)*amf(atom);
                eiveci(atom,dir,phonon,qpt) = eigendi(atom,dir,phonon,qpt)*amf(atom);
                eivec(atom,dir,phonon,qpt) = complex(eivecr(atom,dir,phonon,qpt),eiveci(atom,dir,phonon,qpt));
            end
            eivec_n(atom,:,phonon,qpt) = eivec(atom,:,phonon,qpt)/norm(eivec(atom,:,phonon,qpt));
            for dir = 1:3
                vector(1,(atom-1)*3+dir) = eivec(atom,dir,phonon,qpt);
            end
        end
		for atom = 1:natom
			Erx(atom) = dot(transpose(conj(Rx{atom})),vector);
			Elx(atom) = dot(transpose(conj(Lx{atom})),vector);
			S(phonon,qpt,1,atom) = ((abs(Erx(atom)))^2-(abs(Elx(atom)))^2);
		end		
        chix(phonon,qpt) = sum(S(phonon,qpt,1,:));
		for atom = 1:natom
			Ery(atom) = dot(transpose(conj(Ry{atom})),vector);
			Ely(atom) = dot(transpose(conj(Ly{atom})),vector);
			S(phonon,qpt,2,atom) = ((abs(Ery(atom)))^2-(abs(Ely(atom)))^2);
		end		
        chiy(phonon,qpt) = sum(S(phonon,qpt,2,:));
		for atom = 1:natom
			Erz(atom) = dot(transpose(conj(Rz{atom})),vector);
			Elz(atom) = dot(transpose(conj(Lz{atom})),vector);
			S(phonon,qpt,3,atom) = ((abs(Erz(atom)))^2-(abs(Elz(atom)))^2);
		end		
        chiz(phonon,qpt) = sum(S(phonon,qpt,3,:));
		chit(phonon,qpt) = ((chiz(phonon,qpt))^2 + (chiy(phonon,qpt))^2 + (chix(phonon,qpt))^2)^0.5;
%		clear Sx Sy Sz Ery Erx Erz Ely Elz Elx;
    end
end

%split into contributions from each species

for type = 1:nspecies
    for atom = 1:natom
        if species(atom) == type
            if exist('dummy','var')
                dummy = dummy + S(:,:,:,atom);
            else
                dummy = S(:,:,:,atom);
            end
        end
    end
    Stype{type} = dummy;
    clear dummy;
end

%Calculate gyromagnetic ratios of ions (e/emu)

for atom = 1:natom
    gyro_atom(:,:,atom) = bec(:,:,atom)./(2*mass(species(atom)));
end

%Calculate mode effective charges (e)

mecc = zeros(3,nphonon,nqpt);

for qpt = 1:nqpt
    for phonon = 1:nphonon
        for atom = 1:natom
            mecc(:,phonon,qpt) = mecc(:,phonon,qpt) + bec(:,:,atom)*transpose(eigend(atom,:,phonon,qpt));
            norm(atom,phonon,qpt) = dot(transpose(eigend(atom,:,phonon,qpt)),eigend(atom,:,phonon,qpt));
        end
        mecc_norm(:,phonon,qpt) = mecc(:,phonon,qpt)/sqrt(sum(norm(:,phonon,qpt),1));
        mec_scalar(phonon,qpt) = sum(mecc_norm(:,phonon,qpt));
    end
end

%Calculate magnetic moments of atoms

gyro = zeros(1,3,nphonon,nqpt);

for qpt = 1:nqpt
    for phonon = 1:nphonon
        for atom = 1:natom
           M_atom(phonon,qpt,:,atom) = transpose(gyro_atom(:,:,atom)*(squeeze(S(phonon,qpt,:,atom))))*hbar/nmag;
        end
        M_x(phonon,qpt) = sum(M_atom(phonon,qpt,1,:));
        M_y(phonon,qpt) = sum(M_atom(phonon,qpt,2,:));
        M_z(phonon,qpt) = sum(M_atom(phonon,qpt,3,:));
        M_t(phonon,qpt) = ((M_x(phonon,qpt))^2 + (M_y(phonon,qpt))^2 + (M_z(phonon,qpt))^2)^0.5;  
    end
end

Bfield = M_t*nmag_SI*permeability/(vol/b2A^3*1e-30)*1000;

if plotflag == 1

%make segments through Brillouin zone

nkseg = size(plotlengths,2);
kseg = cell(nkseg,1);
kspecial = zeros(nkseg,1);
nkpts = sum(plotlengths);
nband = size(phf,2);
i = 1;
k = 1;
    
kspecial(k,1) = i;

for j = 1:(plotlengths(1,k)+1)
    kseg{k,1}(j,1)=i;
    i = i + 1;
end
    
for k = 2:nkseg
    kspecial(k,1) = i-1;
    for j = 1:(plotlengths(1,k))
        kseg{k,1}(j,1)=i;
        i = i + 1;
    end
end

kspecial(nkseg+1,1) = nkpts;

%make figure

nbreakval = size(breakval,2);

x = transpose(zeros(nkpts + nbreakval,1));

if nbreakval == 0
    x = [1:nkpts];
else
    p = 0;
    m = 1;
    for n = 1:(nbreakval*2+nkpts)
        if m < (nbreakval + 1)
            if n == (breakval(m) + 2 + (m-1)*2)
                m = m + 1;
            elseif n ~= (breakval(m) + 1 + (m-1)*2)
                p = p + 1;
            end
        else
            p = p +1;
        end
        x(n) = p;
    end
end

clear p m;

z = transpose(zeros(size(x,2),1));

figure
hold on

for j = 1:nband
    sphf = phf(:,j);
    if nbreakval ~= 0
        for m = 1:nbreakval
            sphf = bandbreaker2(sphf,breakval(m)+(m-1)*2);
        end
    end
          
    
     col = real(M_t(j,:)); %Change to plot different chiralities/magnetic moments


    if nbreakval >= 1
        for m = 1:nbreakval
            col = bandbreaker2(col,breakval(m)+(m-1)*2);
        end
    end

    band = surface([x;x],[transpose(sphf);transpose(sphf)],[z;z],[col;col]);
    set(band,'facecol','no','edgecol','interp','LineWidth',2)

    clear col
end

ax = gca;
ax.XTick=kspecial;
ax.XTickLabel=plotlabels;
ax.XGrid='on';
ax.GridAlpha=1;
ax.Box='on';
ax.TickDir='out';
ylabel('Energy / meV'); 
xlim([1,nkpts])

end

function [ x ] = bandbreaker2( x1, breakval )

    x = x1;
    x(breakval+1) = NaN;
    x(breakval+2) = x1(breakval+1);
    for i = (breakval+3):(max(size(x))+2)
        x(i)=x1(i-2);
    end

end


