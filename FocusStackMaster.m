fdir = 'E:\Timelapse-ttr\20220318\20220318';
fname = '20220318_Dex.lif';

mouse = 'Timelapse-ttr';
date = '220211'; %YYMMDD format
run = 1;
ftype = 'sbx_bin';
server = 'Archimedes'; %nickname for server/drive name

bin_t = 3;
bin_xy = 2;

refchannel = 1;
scale = 8;

try
    path = pipe.lab.datapath(mouse,date,run,ftype,server);
    [Nc, Nx, Ny, Nz, Nt, range] = GetDimensions(path,[],[]);
    scale = 4;  
    chunksize = 3; %don't go over 20
    Nchunks = floor(Nt/chunksize);
    anchorref = round(Nchunks/2); %anchor ref must be less than Nchunks
    proj_range = round(0.25*Nz):round(0.75*Nz);
    proj_type = 'mean'; % 'max', 'median', 'mean'
end

path = pipe.lab.datapath(mouse,date,run,ftype,server);

%% if necessary, convert LIF to SBX
if isempty(path)
    reader = bfGetReader(strcat(fdir,filesep,fname));
    M = reader.getMetadataStore();
    Ns = M.getImageCount();
    for i = 1:Ns
        planecount(i) = M.getPlaneCount(i-1);
    end
    [~,SeriesIdx] = max(planecount);
%     SeriesIdx = 3;

    Nc = M.getPixelsSizeC(SeriesIdx - 1).getValue();
    Nx = M.getPixelsSizeX(SeriesIdx - 1).getValue();
    Ny = M.getPixelsSizeY(SeriesIdx - 1).getValue();
    Nz = M.getPixelsSizeZ(SeriesIdx - 1).getValue();
    Nt = M.getPixelsSizeT(SeriesIdx - 1).getValue();

    reader.setSeries(SeriesIdx - 1); %subtract 1 for bioformats 0-indexing

    info = SpoofSBXinfo3D(Ny/bin_xy, Nx/bin_xy, Nz, Nt/bin_t, 1); %generate info file
    savename = strcat(mouse,'_',date,'_',sprintf('%03d',run));
    save(strcat(fdir, filesep, savename,'.mat'),'info', '-v7.3'); %save info mat file
    newpath = strcat(fdir, filesep, savename, '.sbx_bin');
    rw = pipe.io.RegWriter(newpath, info, '.sbx_bin',true);

    w = waitbar(0,'converting to sbx');
    
    nframe = Nz*bin_t;
    
    idx = 1;
    for chunk = 1:bin_t:Nt
        vol = zeros(Nc,Nx,Ny,Nz,bin_t);
        for t = 1:bin_t
            for c = 1:Nc
                for z = 1:Nz
                    if c == 1
                        vol(c,:,:,z,t) = bfGetPlane(reader,idx);
                    end
                    idx = idx+1;
                end
            end
        end
        vol = squeeze(vol(1,:,:,:,:));
        vol = mean(vol,4);
        vol_shrink = imresize3(vol,[Nx/bin_xy,Ny/bin_xy,Nz]);
        rw.write(uint16(vol_shrink));
        waitbar(chunk/(Nt));
    end
    close(w);
end

%% make focused movie
path = pipe.lab.datapath(mouse,date,run,ftype,server);
[Nc, Nx, Ny, Nz, Nt, range] = GetDimensions(path,[],[]);

[fdir, fname, ~] = fileparts(path);
focuspath = strcat(fdir,filesep,fname,'_focus.tif');
entropypath = strcat(fdir,filesep,fname,'_entropy_map.tif');
entmeanpath = strcat(fdir,filesep,fname,'_entropy_heatmap.pdf');

if ~isfile(focuspath)
    ent_width = 10;
    smooth_width = 20;

    %preallocate entropy videos
    mov = zeros(Nx,Ny,Nt);
    entropy_map = zeros(Nx,Ny,Nt);
    entropy_mean = zeros(Nz,Nt);

    H = waitbar(0,'making entropy videos');
    for v_idx = 1:Nt
        vol = pipe.io.sbxRead(path,Nz*(v_idx-1)+1,Nz,1,[]);
        [mov(:,:,v_idx),entropy_map(:,:,v_idx),entropy_mean(:,v_idx)] ...
            = MakeFocusImage(vol,ent_width,smooth_width);
        waitbar(v_idx/Nt);
    end
    close(H);

    %save focused movie
    pipe.io.writeTiff(uint16(mov),focuspath);

    %save entropy movie
    pipe.io.writeTiff(uint16(entropy_map),entropypath);

    %save mean entropy heatmap as heatmap
    figure,imagesc(entropy_mean');
    colormap(turbo);
    saveas(gcf,entmeanpath)
end
%% Load focused movie and high-pass filter
mov_focus = pipe.io.read_tiff(strcat(fdir,filesep,fname,'_focus.tif'));
mov_gauss = imgaussfilt(mov_focus,15);
mov_filt = double(mov_focus) - double(mov_gauss);

%% Do rigid and affine StackReg in Fiji
% stackreg
[tform_rigid, mov_rigid, ~] = StackRegRigid(mov_filt,fdir,'StartSlice',round(Nt/2),'ReturnRegStack',true);
% do affine stackreg
[tform_affine, mov_affine] = StackRegAffine(mov_rigid,fdir,'StartSlice',round(Nt/2),'ReturnRegStack',true);

%% do Demons registration
target = mean(mov_affine(:,:,0.25*Nt:0.75*Nt),3); %make target
w = parfor_progressbar(Nt,'Doing demons registration');
mov_demons = cell(1,Nt);
parfor i = 1:Nt
    slice = mov_affine(:,:,i);
    [D{i},mov_demons{i}] = imregdemons(slice,target,'DisplayWaitbar',false,'AccumulatedFieldSmoothing',3);
    w.iterate(1);
end

mov_demons = [mov_demons{:}];
mov_demons = reshape(mov_demons,size(mov_focus,1),size(mov_focus,2),[]);
implay(rescale(mov_demons));
close(w);
%% Apply warps to raw unfiltered movie

mov_reg = cell(1,Nt);
w = parfor_progressbar(Nt,'applying registrations');
parfor i = 1:Nt
    slice = mov_focus(:,:,i);
    slice_rigid = imwarp(slice,tform_rigid(i),'OutputView',imref2d(size(slice)));
    slice_affine = imwarp(slice_rigid,tform_affine(i),'OutputView',imref2d(size(slice)));
    mov_reg{i} = imwarp(slice_affine,D{i});
%     mov_reg{i} = slice_affine;
    w.iterate(1);
end
mov_reg = [mov_reg{:}];
mov_reg = reshape(mov_reg,size(mov_focus,1),size(mov_focus,2),[]);
close(w);
mov_reg = rescale(mov_reg,0,intmax('uint16'));
implay(rescale(mov_reg));

demonspath = strcat(fdir,filesep,fname,'_demons.tif');
writeTiffFile(uint16(mov_reg),demonspath);

%% Make movie with trace
mov_crop = mov_reg(0.25*Nx:0.75*Ny,0.25*Ny:0.75*Nx,:); %crop edges
trace = squeeze(mean(mean(mov_crop,1),2));

mov_crop = rescale(mov_crop);
ll = prctile(mov_crop(:),1);
ul = prctile(mov_crop(:),99);

mov_crop = rescale(imadjustn(mov_crop,[ll,ul]));
% mov_crop = rescale(mov_crop,ll,ul);
% mov_crop = imadjustn(mov_crop,[0,0.5]);
%%
mov_trace = MakeTraceMov(mov_crop,trace,'Xlabel','Time (hrs)','Ylabel','Fluorescence','FrameInterval',0.5);
v = VideoWriter(strcat(fdir,filesep,fname,'_tracemov.avi'),'Motion JPEG AVI');
v.FrameRate = 15;
open(v);
for i = 1:Nt
    writeVideo(v,mov_trace(i));
end
close(v);