function aseg2wm(aseg)

% Converts aseg.presurf to wm.

wm_labs = [2 41 173 174 175]; % left and right cortical wm, brainstem
gm_labs = [4 11 12 13 26 28 43 50 51 52 58 60]; % sub cortical structures and lat. vent.
% note Natu does not include thalamus labels 10 and 49 in gm_labs, so we also do not


% start with wm as aseg
[p, ~, ~] = fileparts(aseg);
asegd = niftiread(aseg);
wmi = niftiinfo(aseg); % for header info
wm_pre = [p '/wm'];
wm_nii = [wm_pre '.nii'];
wmi.Filename = wm_nii;

wmd = single(zeros(size(asegd)));


for i = wm_labs
    wmd(asegd == i) = 110;
end
clearvars i

for i = gm_labs
    wmd(asegd == i) = 250;
end


niftiwrite(wmd, wm_nii, wmi);

cmd = ['mri_convert -i ' wm_nii ' -o ' wm_pre '.mgz']; % -rt nearest -rl ' aseg
system(cmd);
clearvars cmd
