function ibeat2aseg(ibeat, iFS_dir, out)

% Generates standard FreeSurfer aseg.mgz file from iBEATv2 and Infant FreeSurfer aseg.mgz
% For questions: theodore_turesky@gse.harvard.edu


% build files, filenames, etc.
[~, f, e] = fileparts(ibeat);
fs_ibeat = fullfile(out, ['fs_' f e]);
aseg_in = fullfile(iFS_dir, 'mri', 'aseg.mgz');
aseg_pre = fullfile(out, 'aseg.presurf');
aseg_pre_int = [aseg_pre '_intermed'];
aseg_nii = [aseg_pre '.nii'];
aseg_nii_int = [aseg_pre_int '.nii'];

% ensure same dimensions in ibeat and iFS segmentation files
cmd = ['mri_convert -i ' ibeat ' -o ' fs_ibeat ' -rt nearest -rl ' aseg_in];
system(cmd);
clearvars cmd

% generate nifti version of aseg
cmd = ['mri_convert -i ' aseg_in ' -o ' aseg_nii_int];
system(cmd);
clearvars cmd

% load segmentation (meta)data
ibeat1 = niftiread(fs_ibeat);
asegd = niftiread(aseg_nii_int);
aseg_orig = asegd;
asegi = niftiinfo(aseg_nii_int); % for header info
asegi_int = asegi;
asegi.Filename = aseg_nii;
asegi.Datatype = 'single';
asegi_int.Filename = aseg_nii_int;


% relevant iFS LUT
lw = 2; % cortical wm
rw = 41; % cortical wm
lg = 3; % cortical gm
rg = 42; % cortical gm
lgs = [9 11 12 13 17 18 26 28]; % subcortical gm 9 = thalamus iFS, 10 = thalamus for FS
rgs = [48 50 51 52 53 54 58 60]; % subcortical gm 48 = thalamus for iFS, 49 = thalamus for FS
lwc = 7; % cerebellar wm
rwc = 46; % cerebellar wm
lgc = 8; % cerebellar gm
rgc = 47; % cerebellar gm
lv = 4; % lat ventrical
rv = 43; % lat ventrical
ov = [14 15]; % other ventricals
ot = [172 173 174 175]; % vermis, brainstem

gm = [lg rg lgs rgs lgc rgc ot]; % all GM labels
csf_v = [lv rv ov]; % ventricals and outside brain

% zero skull and outside head
asegd = single(asegd) .* logical(ibeat1);

% generate non-ventrical CSF aseg mask to reduce computational load later
csf_nv = intersect(find(ibeat1 == 1), find(asegd == 0));

% the following are used in the subsequent loop to index ibeat WM (w3), GM (w2),
% and CSF (w1) voxels not assigned in iFS aseg file
w3 = 0;
w2 = 0;
w1 = 0;

% reassign iBEAT voxels with corresponding iFS labels maintaining tissue
% specifications from ibeat
for i = find(ibeat1)'
    switch ibeat1(i)
        case 3
            % label iBEAT WM according to hemisphere
            if ismember(asegd(i), [lw lg lv]) % LH wm, gm, ventricals
                asegd(i) = lw;
            elseif ismember(asegd(i), [rw rg rv]) % RH wm, gm, ventricals
                asegd(i) = rw;
            elseif ismember(asegd(i), [lwc lgc]) % LH cerebellum
                asegd(i) = lwc;  
            elseif ismember(asegd(i), [rwc rgc]) % RH cerebellum
                asegd(i) = rwc;
            elseif ismember(asegd(i), ot) % vermis, brainstem             
                asegd(i) = ot(1, find(ismember(ot, asegd(i))));
            else
                w3 = w3 + 1;
                f_w3(w3,1) = i;                
            end
            
        case 2
            % if iBEAT GM overlaps with GM aseg label, then use aseg label         
            aa = ismember(gm, asegd(i)); 
            if any(aa)
                asegd(i) = gm(1, find(aa));
            else
                w2 = w2 + 1;
                f_w2(w2,1) = i;
            end
            
        case 1
            % label iBEAT CSF according to aseg CSF
            if ismember(asegd(i), lv) % LH ventricals (4)
                asegd(i) = lv;
            elseif ismember(asegd(i), rv) % RH ventricals (43)
                asegd(i) = rv;
            elseif ismember(asegd(i), ov(:,1))
                asegd(i) = ov(:,1);
            elseif ismember(asegd(i), ov(:,2))
                asegd(i) = ov(:,2);
            elseif asegd(i) == 0
                asegd(i) = 0;
            else
                w1 = w1 + 1;
                f_w1(w1,1) = i;
            end

    end
      
     
end

niftiwrite(asegd, aseg_nii_int, asegi_int); % save intermediate aseg file
asegd2 = asegd; % to remove bias from the following process

% generate 3D coordinates for unassigned WM (w3), GM (w2), and CSF (w1) voxels, 
% then identify nearest (Euclidean) label.
if w3
    warning([num2str(w3) ' white matter voxels in iBEAT are not assigned to lateralized WM, GM, or CSF in aseg']);    
    [coord1(:,1), coord1(:,2), coord1(:, 3)] = ind2sub(size(asegd), f_w3); % generate coordinates for unassigned voxels
    coord2 = aseg_labels2coords(asegd, [lw rw lwc rwc ot]); % generate coordinates for relevant  aseg labels
    [k, dist] = dsearchn(coord2(:, 1:3), coord1); % label according to nearest WM label  
    asegd2(f_w3) = coord2(k, 4);
    clearvars coord1 coord2 k dist
end

if w2
    warning([num2str(w2) ' gray matter voxels in iBEAT are not assigned to aseg GM']);
    [coord1(:,1), coord1(:,2), coord1(:,3)] = ind2sub(size(asegd), f_w2); % generate coordinates for unassigned voxels
    coord2 = aseg_labels2coords(asegd, gm); % generate coordinates for relevant aseg labels
    [k, dist] = dsearchn(coord2(:,1:3), coord1); % label according to nearest GM label    
    asegd2(f_w2) = coord2(k, 4);
    clearvars coord1 coord2 k dist
end

if w1
    warning([num2str(w1) ' CSF voxels in iBEAT are not assigned to aseg CSF']);
     [coord1(:,1), coord1(:,2), coord1(:,3)] = ind2sub(size(asegd), f_w1); % generate coordinates for unassigned voxels
     coord2v = aseg_labels2coords(asegd, csf_v); % generate coordinates for relevant  aseg labels
     [coord2nv(:,1), coord2nv(:,2), coord2nv(:,3)] = ind2sub(size(asegd), csf_nv);
     coord2 = vertcat(coord2v, [coord2nv zeros(size(coord2nv,1),1)]);
     [k, dist] = dsearchn(coord2(:,1:3), coord1); % label according to nearest CSF label   
     asegd2(f_w1) = coord2(k, 4);
     clearvars coord1 coord2 k dist
end

% Because iFS seemed to perform better on subcortical areas, relabel
% all subcortical (non-cerebellum areas) with aseg
for ii = [lgs rgs]
    asegd2(aseg_orig == ii) = ii;
end

% subcortical gm 9 = thalamus iFS, 10 = thalamus for FS
% subcortical gm 48 = thalamus for iFS, 49 = thalamus for FS

asegd2(asegd2 == 9) = 10;
asegd2(asegd2 == 48) = 49;

niftiwrite(asegd2, aseg_nii, asegi);

% save ibeat-adjusted aseg file
cmd = ['mri_convert -i ' aseg_nii ' -o ' aseg_pre '.mgz']; % -rt nearest -rl ' aseg
system(cmd);
clearvars cmd
