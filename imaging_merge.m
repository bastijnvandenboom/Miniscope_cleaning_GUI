
% changed by bastijn number of times
% 1) changed figure: easier to use now
% 2) when you press Add Lazy, you go to the next pair
% 3) when you press Add to delete list, you go to the next neuron
% 4) show number of cells after merging5


function varargout = imaging_merge(varargin)
% IMAGING_MERGE MATLAB code for imaging_merge.fig
%      IMAGING_MERGE, by itself, creates a new IMAGING_MERGE or raises the existing
%      singleton*.
%
%      H = IMAGING_MERGE returns the handle to a new IMAGING_MERGE or the handle to
%      the existing singleton*.
%
%      IMAGING_MERGE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGING_MERGE.M with the give0n input arguments.
%
%      IMAGING_MERGE('Property','Value',...) creates a new IMAGING_MERGE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before imaging_merge_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to imaging_merge_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help imaging_merge

% Last Modified by GUIDE v2.5 18-Sep-2019 15:57:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @imaging_merge_OpeningFcn, ...
    'gui_OutputFcn',  @imaging_merge_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before imaging_merge is made visible.
function imaging_merge_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to imaging_merge (see VARARGIN)

% Choose default command line output for imaging_merge
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes imaging_merge wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = imaging_merge_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in LoadFile.
function LoadFile_Callback(hObject, eventdata, handles)
% hObject    handle to LoadFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename1,filepath1]=uigetfile({'*.*','All Files'},...
    'Select Data File 1');
handles.rawdata1=load([filepath1 filename1]);
handles.N_neurons = size(handles.rawdata1.results.C_raw,1);
handles.updatedresults = handles.rawdata1.results.C_raw;
handles.updatedC = handles.rawdata1.results.C;
handles.updatedA = handles.rawdata1.results.A;
handles.mergecells = [];
handles.delcells=[];
set(handles.cell1,'String',[]);
set(handles.cell2,'String',[]);
set(handles.deletelist,'String',[]);
guidata(hObject,handles)




% --- Executes on button press in pre.
function pre_Callback(hObject, eventdata, handles)
% hObject    handle to pre (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
C = handles.updatedresults;
if handles.currentcompare~=1
    k = handles.currentcompare-1;
    handles.currentcompare = k;
    plot(handles.ax1,C(handles.r(k),:),'-r');
    hold(handles.ax1,'on')
    plot(handles.ax1,C(handles.co(k),:),'-b');
    hold(handles.ax1,'off')
    str = sprintf('Cell 1 - %d \n Cell 2 - %d \n Dist b/w them - %0.2d Crosscorr - %0.2d', handles.r(k),handles.co(k),handles.pwdist(handles.r(k),handles.co(k)),handles.crosscoef(handles.r(k),handles.co(k)));
    title(handles.ax1,str);
    legend(handles.ax1,'Cell 1','Cell 2');
    cc = tril(handles.crosscoef);
    pw = reshape(handles.pwdist,1,size(handles.pwdist,1)*size(handles.pwdist,2));
    cc = reshape(cc,1,size(cc,1)*size(cc,2));
    ind0 = find(pw==0);
    pw(ind0)=[];
    cc(ind0)=[];
    axes(handles.ax2);
    imshow(handles.rawdata1.results.Cn,[0 1]);
    hold(handles.ax2,'on')
    
    img = handles.rawdata1.results.Cn;
    ind = [handles.r(handles.currentcompare) handles.co(handles.currentcompare)];
    thr =  str2num(get(handles.thr,'String'));
    with_label = false;
    Coor = get_contours(handles,thr,ind,1);
    handles.Coor = Coor;
    d1 = handles.rawdata1.results.options.d1;
    d2 = handles.rawdata1.results.options.d2;
    plot_contours(handles.updatedA(:, ind),d1,d2, img, thr,with_label, [], handles.Coor, 2);
    colormap gray;
    hold(handles.ax2,'off')
    % x_pixel = handles.rawdata1.results.options.d1;
    % y_pixel = handles.rawdata1.results.options.d2;
    % xlim(handles.ax2,[0 x_pixel])
    % ylim(handles.ax2,[0 y_pixel])
    title(handles.ax2,'Spatial locations of the cells')
    plot(handles.distcorr,pw,cc,'.r','LineStyle','none','MarkerSize',10);
    hold(handles.distcorr,'on')
    plot(handles.distcorr,handles.pwdist(handles.r(handles.currentcompare),handles.co(handles.currentcompare)),handles.crosscoef(handles.r(handles.currentcompare),handles.co(handles.currentcompare)),'.k','LineStyle','none','MarkerSize',15)
    hold(handles.distcorr,'off')
end
guidata(hObject,handles)




% --- Executes on button press in next.
function next_Callback(hObject, eventdata, handles)
% hObject    handle to next (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
C = handles.rawdata1.results.C_raw;

if handles.currentcompare~=length(handles.r)
    k = handles.currentcompare+1;
    handles.currentcompare = k;
    plot(handles.ax1,C(handles.r(k),:),'-r');
    hold(handles.ax1,'on')
    plot(handles.ax1,C(handles.co(k),:),'-b');
    hold(handles.ax1,'off')
    str = sprintf('Cell 1 - %d \n Cell 2 - %d \n Dist b/w them - %0.2d Crosscorr - %0.2d', handles.r(k),handles.co(k),handles.pwdist(handles.r(k),handles.co(k)),handles.crosscoef(handles.r(k),handles.co(k)));
    title(handles.ax1,str);
    legend(handles.ax1,'Cell 1','Cell 2');
    cc = tril(handles.crosscoef);
    pw = reshape(handles.pwdist,1,size(handles.pwdist,1)*size(handles.pwdist,2));
    cc = reshape(cc,1,size(cc,1)*size(cc,2));
    ind0 = find(pw==0);
    pw(ind0)=[];
    cc(ind0)=[];
    axes(handles.ax2);
    imshow(handles.rawdata1.results.Cn,[0 1]);
    
    img = handles.rawdata1.results.Cn;
    ind = [handles.r(handles.currentcompare) handles.co(handles.currentcompare)];
    thr =  str2num(get(handles.thr,'String'));
    with_label = false;
    Coor = get_contours(handles,thr,ind,1);
    handles.Coor = Coor;
    d1 = handles.rawdata1.results.options.d1;
    d2 = handles.rawdata1.results.options.d2;
    plot_contours(handles.updatedA(:, ind),d1,d2, img, thr,with_label, [], handles.Coor, 2);
    colormap gray;
    hold(handles.ax2,'off')
    % x_pixel = handles.rawdata1.results.options.d1;
    % y_pixel = handles.rawdata1.results.options.d2;
    % xlim(handles.ax2,[0 x_pixel])
    % ylim(handles.ax2,[0 y_pixel])
    title(handles.ax2,'Spatial locations of the cells')
    plot(handles.distcorr,pw,cc,'.r','LineStyle','none','MarkerSize',10);
    hold(handles.distcorr,'on')
    plot(handles.distcorr,handles.pwdist(handles.r(handles.currentcompare),handles.co(handles.currentcompare)),handles.crosscoef(handles.r(handles.currentcompare),handles.co(handles.currentcompare)),'.k','LineStyle','none','MarkerSize',15)
    hold(handles.distcorr,'off')
end
guidata(hObject,handles)


% --- Executes on selection change in cell1.
function cell1_Callback(hObject, eventdata, handles)
% hObject    handle to cell1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns cell1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cell1

set(handles.cell1,'String',unique([handles.r handles.co]));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function cell1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cell1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in cell2.
function cell2_Callback(hObject, eventdata, handles)
% hObject    handle to cell2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns cell2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cell2
set(handles.cell2,'String',unique([handles.r handles.co]));
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function cell2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cell2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in merge.
function merge_Callback(hObject, eventdata, handles)
% hObject    handle to merge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
k1 = get(handles.mergelist,'String');
tmp1=[];
tmp2=[];
for i = 1:length(k1)
    tmp3 = sort(str2num(cell2mat(k1(i))),'ascend');
    tmp4=[];
for u = 1:length(tmp3)
    tmp4 = [tmp4;handles.updatedresults(tmp3(u),:)];
end
    merged_signal(i,:) = mean(tmp4);
    merged_contour(:,i) = handles.updatedA(:,tmp3(1));
    tmp1 = [tmp1 tmp3(1)];
    tmp2 = [tmp2 tmp3(2:end)];
end
handles.updatedresults(tmp1,:) = merged_signal;
handles.updatedresults(tmp2,:)=[];
handles.updatedA(:,tmp1) = merged_contour;
handles.updatedA(:,tmp2)=[];
set(handles.mergelist,'String',[]);
% set(handles.mergecells,'String',[]);
set(handles.merge_alert,'String','Done merging! Restart the gui if you need to redo merging')
handles.pwdist=[];
handles.crosscoef=[];
handles.r=[];
handles.co=[];
handles.mergecells=[];

% show number of cells
str = sprintf('Updated number of neurons - %d',size(handles.updatedresults,1));
set(handles.N_final,'String',str);

guidata(hObject,handles);



% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

folder_name = uigetdir;
file_str = sprintf('%s\\updated_Craw.mat',folder_name);
results = handles.rawdata1.results;
results.C_raw = handles.updatedresults;
results.A =  handles.updatedA;
fig1_str = sprintf('%s\\corrVsdist.png',folder_name);

saveas(handles.distcorr,fig1_str);

save(file_str,'results');
save_str = sprintf('Saved at %s',file_str);
set(handles.save_text,'String',save_str);
axes(handles.ax2);
imshow(handles.rawdata1.results.Cn,[0 1]);

img = handles.rawdata1.results.Cn;
ind = 1:size(handles.updatedA,2);
thr =  str2num(get(handles.thr,'String'));
with_label = false;
Coor = get_contours(handles,thr,ind,1);
handles.Coor = Coor;
d1 = handles.rawdata1.results.options.d1;
d2 = handles.rawdata1.results.options.d2;
plot_contours(handles.updatedA(:, ind),d1,d2, img, thr,with_label, [], handles.Coor, 2);
colormap gray;
hold(handles.ax2,'off')
guidata(hObject,handles);






function dist_Callback(hObject, eventdata, handles)
% hObject    handle to dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dist as text
%        str2double(get(hObject,'String')) returns contents of dist as a double
handles.distthr =  str2double(get(handles.dist,'String'))


% --- Executes during object creation, after setting all properties.
function dist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function corr_Callback(hObject, eventdata, handles)
% hObject    handle to corr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of corr as text
%        str2double(get(hObject,'String')) returns contents of corr as a double
handles.corrthr =  str2double(get(handles.corr,'String'))


% --- Executes during object creation, after setting all properties.
function corr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to corr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function thr_Callback(hObject, eventdata, handles)
% hObject    handle to thr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thr as text
%        str2double(get(hObject,'String')) returns contents of thr as a double


% --- Executes during object creation, after setting all properties.
function thr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function snr_Callback(hObject, eventdata, handles)
% hObject    handle to snr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of snr as text
%        str2double(get(hObject,'String')) returns contents of snr as a double


% --- Executes during object creation, after setting all properties.
function snr_CreateFcn(hObject, eventdata, handles)
% hObject    handle to snr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in deletelist.
function deletelist_Callback(hObject, eventdata, handles)
% hObject    handle to deletelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns deletelist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from deletelist
k1 = get(handles.deletelist,'Value');
k2 = str2num(get(handles.deletelist,'String'));
tmp = k2(k1);

C=  handles.updatedresults;
plot(handles.ax1,C(tmp,:));


imshow(handles.rawdata1.results.Cn,[0 1]);
hold(handles.ax2,'on')

img = handles.rawdata1.results.Cn;
ind = tmp;
thr =  str2num(get(handles.thr,'String'));
with_label = false;
Coor = get_contours(handles,thr,ind,1);
handles.Coor = Coor;
d1 = handles.rawdata1.results.options.d1;
d2 = handles.rawdata1.results.options.d2;
plot_contours(handles.updatedA(:, ind),d1,d2, img, thr,with_label, [], handles.Coor, 2);
colormap gray;
hold(handles.ax2,'off')
% x_pixel =handles.rawdata1.results.options.d1;
% y_pixel =handles.rawdata1.results.options.d2;
% xlim(handles.ax2,[0 x_pixel])
% ylim(handles.ax2,[0 y_pixel])
title(handles.ax2,'Spatial locations of the cells');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function deletelist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to deletelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in del.
function del_Callback(hObject, eventdata, handles)
% hObject    handle to del (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
k1 = str2num(get(handles.deletelist,'String'));
handles.updatedresults(k1,:)=[];
handles.updatedA(:,k1)=[];
str = sprintf('Updated number of neurons - %d',size(handles.updatedresults,1));
set(handles.N_final,'String',str);
handles.pwdist = [];
handles.crosscoef=[];
handles.r=[];
handles.co=[];
guidata(hObject,handles);





% --- Executes on selection change in celllist.
function celllist_Callback(hObject, eventdata, handles)
% hObject    handle to celllist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns celllist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from celllist
set(handles.celllist,'String',1:size(handles.updatedresults,1));
k= get(handles.celllist,'Value');
tmp = handles.updatedresults;
plot(handles.ax1,tmp(k,:),'-r');
str = sprintf(' SNR- %0.2d',handles.SNR(k));
title(handles.ax1,str);
axes(handles.ax2)
imshow(handles.rawdata1.results.Cn,[0 1]);
hold(handles.ax2,'on')

img = handles.rawdata1.results.Cn;
ind = k;
thr =  str2num(get(handles.thr,'String'));
with_label = false;
Coor = get_contours(handles,thr,ind,1);
handles.Coor = Coor;
d1 = handles.rawdata1.results.options.d1;
d2 = handles.rawdata1.results.options.d2;
plot_contours(handles.updatedA(:, ind),d1,d2, img, thr,with_label, [], handles.Coor, 2);
colormap gray;
hold(handles.ax2,'off')
% x_pixel = handles.rawdata1.results.options.d1;
% y_pixel = handles.rawdata1.results.options.d2;
% xlim(handles.ax2,[0 x_pixel])
% ylim(handles.ax2,[0 y_pixel])
title(handles.ax2,'Spatial locations of the cells');







% --- Executes during object creation, after setting all properties.
function celllist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to celllist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in dellist.
function dellist_Callback(hObject, eventdata, handles)
% hObject    handle to dellist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
k1 = handles.k;
tmp_list = handles.delcells;
tmp_list = [tmp_list k1];
handles.delcells = tmp_list;
set(handles.deletelist,'String',tmp_list,'Value',length(tmp_list));

% add to go to next neuron (unless it's the last neuron)
if handles.k ~= size(handles.updatedresults,1)
handles.k= handles.k + 1;
set(handles.slidertoaddtodellist,'Max',size(handles.updatedresults,1),'Value',handles.k,'SliderStep',[1/(size(handles.updatedresults,1)-1),1]);

tmp = handles.updatedresults;
plot(handles.ax1,tmp(handles.k,:),'-r');
str = sprintf(' SNR- %0.2d Cell no - %d',handles.SNR(handles.k),handles.k);
title(handles.ax1,str);
axes(handles.ax2)
imshow(handles.rawdata1.results.Cn,[0 1]);
hold(handles.ax2,'on')

img = handles.rawdata1.results.Cn;
ind = handles.k;
thr =  str2num(get(handles.thr,'String'));
with_label = false;
Coor = get_contours(handles,thr,ind,1);
handles.Coor = Coor;
d1 = handles.rawdata1.results.options.d1;
d2 = handles.rawdata1.results.options.d2;
plot_contours(handles.updatedA(:, ind),d1,d2, img, thr,with_label, [], handles.Coor, 2);
colormap gray;
hold(handles.ax2,'off')
% x_pixel =handles.rawdata1.results.options.d1;
% y_pixel = handles.rawdata1.results.options.d2;
% xlim(handles.ax2,[0 x_pixel])
% ylim(handles.ax2,[0 y_pixel])
title(handles.ax2,'Spatial locations of the cells');
end

guidata(hObject,handles);




function [CC,jsf] = plot_contours(Aor,d1,d2,Cn,thr,display_numbers,max_number,Coor, ln_wd, ln_col)

% save and plot the contour traces of the found spatial components againsts
% specified background image. The contour is drawn around the value above
% which a specified fraction of energy is explained (default 99%)

if nargin < 5 || isempty(max_number)
    max_number = size(Aor,2);
else
    max_number = min(max_number,size(Aor,2));
end
if nargin < 4 || isempty(display_numbers)
    display_numbers = 0;
end
if nargin < 3 || isempty(thr)
    thr = 0.995;
end
if ~exist('ln_wd', 'var') || isempty(ln_wd)
    ln_wd = 1; % linewidth;
end
units = 'centimeters';
fontname = 'helvetica';

%fig3 = figure;
%     set(gcf, 'PaperUnits', units,'Units', units)
%     set(gcf, 'PaperPosition',[5, 5, 12, 12])
%     set(gcf, 'Position',3*[5, 5, 12, 12])

%set(gca,'XTick',[],'YTick',[]);
posA = get(gca,'position');
set(gca,'position',posA);
%cbar = colorbar('south','TickDirection','out');
if (0)
    cbar = colorbar('TickDirection','out');
    cpos = get(cbar,'position');
    %cpos = [posA(1),posA(2)-cpos(4)-0.01,posA(3),cpos(4)];
    ylabel(cbar,'Average neighbor correlation');
    set(cbar,'position',cpos,'TickDirection','in');
    set(cbar,'fontweight','bold','fontsize',14,'fontname',fontname);
end
%hold on; scatter(cm(:,2),cm(:,1),'ko'); hold off;
%v = axis;
%handle = title('Correlation image and identified spatial footprints','fontweight','bold','fontsize',14,'fontname',fontname);
hold on;
if ~exist('cmap', 'var') || isempty(cmap)
    cmap = hot(3*size(Aor,2));
else
    cmap = repmat(reshape(ln_col, 1, []), size(Aor,2), 1);
end
if ~(nargin < 6 || isempty(Coor))
    CC = Coor;
    for i = 1:size(Aor,2)
        %         cont = medfilt1(Coor{i}')';
        cont = Coor{i};
        if size(cont,2) > 1
            plot(cont(1,1:end),cont(2,1:end),'Color',cmap(i+size(Aor,2),:), 'linewidth', ln_wd); hold on;
        end
    end
else
    CC = cell(size(Aor,2),1);
    CR = cell(size(Aor,2),2);
    for i = 1:size(Aor,2)
        A_temp = full(reshape(Aor(:,i),d1,d2));
        A_temp = medfilt2(A_temp,[3,3]);
        A_temp = A_temp(:);
        [temp,ind] = sort(A_temp(:).^2,'ascend');
        temp =  cumsum(temp);
        ff = find(temp > (1-thr)*temp(end),1,'first');
        if ~isempty(ff)
            CC{i} = contour(reshape(A_temp,d1,d2),[0,0]+A_temp(ind(ff)),'LineColor',cmap(i+size(Aor,2),:), 'linewidth', ln_wd);
            fp = find(A_temp >= A_temp(ind(ff)));
            [ii,jj] = ind2sub([d1,d2],fp);
            CR{i,1} = [ii,jj]';
            CR{i,2} = A_temp(fp)';
        end
        hold on;
    end
end

if display_numbers
    cm = com(Aor(:,1:end),d1,d2);
    lbl = strtrim(cellstr(num2str((1:size(Aor,2))')));
    text((cm(1:max_number,2)),(cm(1:max_number,1)),lbl(1:max_number),'color',[0,0,0],'fontsize',16,'fontname',fontname,'fontweight','bold');
end
axis off;
if ~(nargin < 6 || isempty(Coor))
    jsf = [];
else
    for i = 1:size(Aor,2);
        if ~isempty(CR{i,1})
            jsf(i) = struct('id',i,...
                'coordinates',CR{i,1}',...
                'values',CR{i,2},...
                'bbox',[min(CR{i,1}(1,:)),max(CR{i,1}(1,:)),min(CR{i,1}(2,:)),max(CR{i,1}(2,:))],...
                'centroid',cm(i,:));
        end
        if i == 1
            jsf = repmat(jsf,size(Aor,2),1);
        end
    end
end





function [b, sn] = estimate_baseline_noise(y, bmin)
%% estimate the baseline and noise level for single calcium trace
%% inputs:
%   y:  T*1 vector, calcium trace
%   thr: scalar, threshold for choosing points for fitting a gaussian
%   distribution.
%   bmin: scalar, minimum value of the baseline, default(-inf)
%% outputs:
%   b: scalar, baseline
%   sn: scalar, sigma of the noise
%% Author: Pengcheng Zhou, Carnegie Mellon University, 2016

%% input arguments
if ~exist('bmin', 'var') || isempty(bmin)
    bmin = -inf;
end

%% create the histogram for fitting
temp = quantile(y, 0:0.1:1);
dbin = max(min(diff(temp))/3, (max(temp)-min(temp))/1000);
bins = temp(1):dbin:temp(end);
nums = hist(y, bins);

%% fit a gaussian distribution: nums = A * exp(-(bins-b)/(2*sig^2))
if isempty(bins)
    b = mean(y);
    sn = 0;
    return;
end
[b, sn] = fit_gauss1(bins, nums, 0.3, 3);

if b<bmin
    b = bmin;
    sn = fit_gauss1(bins-bmin, nums, 0.3, 3,false );
end

function Coor = get_contours(h2, thr, ind,updated)
   if updated
        A_ = h2.updatedA;
        A_ = A_(:,ind);
    else
          A_ = h2.rawdata1.results.A;

              A_ = A_(:, ind);
    end


          num_neuron = size(A_,2);
          if num_neuron==0
              Coor ={};
              return;
          else
              Coor = cell(num_neuron,1);
          end
          d1 = h2.rawdata1.results.options.d1,1;
          d2 = h2.rawdata1.results.options.d2,2;
          %             tmp_kernel = strel('square', 3);
          for m=1:num_neuron
              % smooth the image with median filter
              A_temp = reshape(full(A_(:, m)),d1,d2);
              % find the threshold for detecting nonzero pixels

              A_temp = A_temp(:);
              [temp,ind] = sort(A_temp(:).^2,'ascend');
              temp =  cumsum(temp);
              ff = find(temp > (1-thr)*temp(end),1,'first');
              thr_a = A_temp(ind(ff));
              A_temp = reshape(A_temp,d1,d2);

              % crop a small region for computing contours
              [tmp1, tmp2, ~] = find(A_temp);
              if isempty(tmp1)
                  Coor{m} = zeros(2,1);
                  continue;
              end
              rmin = max(1, min(tmp1)-3);
              rmax = min(d1, max(tmp1)+3);
              cmin = max(1, min(tmp2)-3);
              cmax = min(d2, max(tmp2)+3);
              A_temp = A_temp(rmin:rmax, cmin:cmax);

              if nnz(A_temp)>36
                  l = bwlabel(medfilt2(A_temp>thr_a));
              else
                  l = bwlabel(A_temp>=thr_a);
              end
              l_most = mode(l(l>0));
              if isnan(l_most)
                  Coor{m} = zeros(2, 1);
                  continue;
              end
              ind = (l==l_most);
              A_temp(ind) =  max(A_temp(ind), thr_a);
              A_temp(~ind) = min(A_temp(~ind), thr_a*0.99);

              pvpairs = { 'LevelList' , thr_a, 'ZData', A_temp};
              h = matlab.graphics.chart.primitive.Contour(pvpairs{:});
              temp = h.ContourMatrix;
              if isempty(temp)
                  temp = get_contours(hObject,(thr+1)/2, ind_show(m));
                  Coor{m} = temp{1};
                  continue;
              else
                  temp(:, 1) = temp(:, 2);
                  temp = medfilt1(temp')';
                  temp(:, 1) = temp(:, end);
                  Coor{m} = bsxfun(@plus, temp, [cmin-1; rmin-1]);
              end
end



function [mu, sig, A] = fit_gauss1(x, y, thr, maxIter, mu_fix)
%% fit a gaussin curve given the data (x, y): y = A*exp(-(x-mu)^2/2/sig^2)

%% inputs:
%   x: T*1 vector
%   y: T*1 vector
%   thr: scalar, threshold for selecting points for fitting
%   maxIter: scalar, maximum iteration
%   mu_fix: boolean, fix the center to be 0 or not
%% outputs:
%   mu: scalar, center of the peak
%   sig: scalar, gaussian width
%   A:  scalar, amplitude of the gaussian curve

%% Author: Pengcheng Zhou, Carnegie Mellon University, 2016
%% Reference:
%   A Simple Algorithm for Fitting a Gaussian Function , Hongwei Guo, 2011

%% input arguments
x = reshape(x, [], 1);
y = reshape(y, [], 1);
T = length(y);

if ~exist('thr', 'var') || isempty(thr)
    thr = 0.1;
end
if ~exist('maxIter', 'var') || isempty(maxIter)
    maxIter = 5;
end
if ~exist('mu_fix', 'var') || isempty(mu_fix)
    mu_fix = false;
end
ind = (y>max(y)*thr);
x = x(ind);
y = y(ind);

x2 = x.^2;
x3 = x.^3;
x4 = x.^4;
y2 = y.^2;
logy = log(y);
y2logy = y2.*logy;
vec1 = ones(1, length(y));

warning('off','MATLAB:nearlySingularMatrix');
warning('off','MATLAB:SingularMatrix');
%% fit the curve
if mu_fix % fix the mu to be 0
    for miter=1:maxIter
        M = [vec1*y2, x2'*y2; ...
            x2'*y2, x4'*y2];
        b = [vec1*y2logy; x2'*y2logy];
        p = M\b;
        
        logy = p(1)*vec1' + p(2)*x2;
        y = exp(logy);
        y2 = y.^2;
        y2logy = y2.*logy;
    end
    
    mu= 0;
    sig = sqrt(-0.5/p(2));
    A = exp(p(1));
else
    for miter=1:maxIter
        M = [vec1*y2, x'*y2, x2'*y2; ...
            x'*y2, x2'*y2, x3'*y2; ...
            x2'*y2, x3'*y2, x4'*y2];
        b = [vec1*y2logy; x'*y2logy; x2'*y2logy];
        p = (M)\b;
        
        logy = p(1)*vec1' + p(2)*x + p(3)*x2;
        y = exp(logy);
        y2 = y.^2;
        y2logy = y2.*logy;
    end
    mu= -p(2)/2/p(3);
    sig = abs(sqrt(-0.5/p(3)));
    A = exp(p(1)-0.25*p(2)^2/p(3));
end

warning('on','MATLAB:nearlySingularMatrix');
warning('on','MATLAB:SingularMatrix');


% --- Executes on slider movement.
function slidertoaddtodellist_Callback(hObject, eventdata, handles)
% hObject    handle to slidertoaddtodellist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(handles.slidertoaddtodellist,'Max',size(handles.updatedresults,1),'Value',round(get(handles.slidertoaddtodellist,'Value')),'SliderStep',[1/(size(handles.updatedresults,1)-1),1]);
handles.k= get(handles.slidertoaddtodellist,'Value');
tmp = handles.updatedresults;
plot(handles.ax1,tmp(handles.k,:),'-r');
str = sprintf(' SNR- %0.2d Cell no - %d',handles.SNR(handles.k),handles.k);
title(handles.ax1,str);
axes(handles.ax2)
imshow(handles.rawdata1.results.Cn,[0 1]);
hold(handles.ax2,'on')

img = handles.rawdata1.results.Cn;
ind = handles.k;
thr =  str2num(get(handles.thr,'String'));
with_label = false;
Coor = get_contours(handles,thr,ind,1);
handles.Coor = Coor;
d1 = handles.rawdata1.results.options.d1;
d2 = handles.rawdata1.results.options.d2;
plot_contours(handles.updatedA(:, ind),d1,d2, img, thr,with_label, [], handles.Coor, 2);
colormap gray;
hold(handles.ax2,'off')
% x_pixel =handles.rawdata1.results.options.d1;
% y_pixel = handles.rawdata1.results.options.d2;
% xlim(handles.ax2,[0 x_pixel])
% ylim(handles.ax2,[0 y_pixel])
title(handles.ax2,'Spatial locations of the cells');

guidata(hObject,handles);




% --- Executes during object creation, after setting all properties.
function slidertoaddtodellist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slidertoaddtodellist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



% --- Executes on selection change in mergelist.
function mergelist_Callback(hObject, eventdata, handles)
% hObject    handle to mergelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns mergelist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from mergelist
k1 = get(handles.mergelist,'Value');
k2 = get(handles.mergelist,'String');
tmp = k2(k1);
tmp = str2num(cell2mat(tmp));
C=  handles.rawdata1.results.C_raw;
for u = 1:length(tmp)
plot(handles.ax1,C(tmp(u),:));
hold(handles.ax1,'on')
end
hold(handles.ax1,'off')

imshow(handles.rawdata1.results.Cn,[0 1]);
hold(handles.ax2,'on')

img = handles.rawdata1.results.Cn;
ind = tmp;
thr =  str2num(get(handles.thr,'String'));
with_label = false;
Coor = get_contours(handles,thr,ind,1);
handles.Coor = Coor;
d1 = handles.rawdata1.results.options.d1;
d2 = handles.rawdata1.results.options.d2;
plot_contours(handles.updatedA(:, ind),d1,d2, img, thr,with_label, [], handles.Coor, 2);
colormap gray;
hold(handles.ax2,'off')
% x_pixel =handles.rawdata1.results.options.d1;
% y_pixel =handles.rawdata1.results.options.d2;
% xlim(handles.ax2,[0 x_pixel])
% ylim(handles.ax2,[0 y_pixel])
title(handles.ax2,'Spatial locations of the cells');
guidata(hObject,handles);




% --- Executes during object creation, after setting all properties.
function mergelist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mergelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in add2merge.
function add2merge_Callback(hObject, eventdata, handles)
% hObject    handle to add2merge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
C = handles.rawdata1.results.C_raw;
un_cells = unique([handles.r handles.co]);
cell_1 = un_cells(cell2mat(get(handles.cell1,'Value')));
cell_2 = un_cells(cell2mat(get(handles.cell2,'Value')));
m = handles.mergecells;
m = [m; {mat2str([cell_1 cell_2])} ];
set(handles.mergelist,'String',m,'Value',length(m))
handles.mergecells=m;
guidata(hObject,handles)




% --- Executes on button press in delete4mmerge.
function delete4mmerge_Callback(hObject, eventdata, handles)
% hObject    handle to delete4mmerge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
k1 =(get(handles.mergelist,'Value'));
k2 =(get(handles.mergelist,'String'));
k2(k1,:)=[];
handles.mergecells(k1,:)=[];
set(handles.mergelist,'String',k2,'Value',size(k2,1));
guidata(hObject,handles)


% --- Executes on button press in rem4mdellist.
function rem4mdellist_Callback(hObject, eventdata, handles)
% hObject    handle to rem4mdellist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
k1 =(get(handles.deletelist,'Value'));
k2 =(get(handles.deletelist,'String'));
k2(k1,:)=[];
handles.delcells(k1)=[];
set(handles.deletelist,'String',k2,'Value',size(k2,1));
guidata(hObject,handles)




% --- Executes on button press in update.
function update_Callback(hObject, eventdata, handles)
% hObject    handle to update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.SNR=[];
for i = 1:size(handles.rawdata1.results.C_raw,1)
    [pk,loc] = findpeaks(handles.updatedC(i,:),3);
     P{i} = loc;
    Pks{i}= pk;
end
for i=1:size(handles.updatedresults,1)
    
    [b, sn] = estimate_baseline_noise(handles.updatedresults(i,:));
     handles.SNR(i) = median(Pks{i})/sn;
%     handles.SNR(i) = Pks(i)/sn;
    clear sn;
end
ind_snr = find(handles.SNR < str2num(get(handles.snr,'String')));
handles.updatedresults(ind_snr,:)=[];
handles.updatedA(:,ind_snr)=[];
 handles.N_neurons = size(handles.updatedresults,1);
 A = reshape(full(handles.updatedA)',handles.N_neurons,handles.rawdata1.results.options.d1,handles.rawdata1.results.options.d2);
% handles.N_neurons = size(handles.updatedresults,1);
% A(ind_snr,:,:) = [];
% handles.updatedA = reshape(A,handles.rawdata1.results.options.d1*handles.rawdata1.results.options.d2,handles.N_neurons);
handles.SNR(ind_snr)=[];
Cent_N = zeros(size(handles.updatedresults,1),2);


for i = 1:handles.N_neurons
    [tmp_r,tmp_c] = find(squeeze(A(i,:,:)));
    Cent_N(i,:) =[mean(tmp_r) mean(tmp_c)];
end

C = handles.updatedresults;
for i = 1:handles.N_neurons
    for j = 1:handles.N_neurons
        handles.pwdist(i,j) = sqrt((Cent_N(i,1) - Cent_N(j,1))^2 + (Cent_N(i,2) - Cent_N(j,2))^2);
        handles.crosscoef(i,j) = corr(C(i,:)',C(j,:)');
    end
end
handles.pwdist = tril(handles.pwdist);
cc = tril(handles.crosscoef);
pw = reshape(handles.pwdist,1,size(handles.pwdist,1)*size(handles.pwdist,2));
cc = reshape(cc,1,size(cc,1)*size(cc,2));
ind0 = find(pw==0);
pw(ind0)=[];
cc(ind0)=[];

[handles.r,handles.co] = find(handles.pwdist>0 & handles.pwdist < str2num(get(handles.dist,'String')) & handles.crosscoef>str2num(get(handles.corr,'String')));
plot(handles.ax1,C(handles.r(1),:),'-r');
hold(handles.ax1,'on')
plot(handles.ax1,C(handles.co(1),:),'-b');
hold(handles.ax1,'off')
set(handles.cell1,'String',handles.r);
set(handles.cell2,'String',handles.co)
str = sprintf('Cell 1 - %d \n Cell 2 - %d \n Dist b/w them - %0.2d Crosscorr - %0.2d', handles.r(1),handles.co(1),handles.pwdist(handles.r(1),handles.co(1)),handles.crosscoef(handles.r(1),handles.co(1)));
title(handles.ax1,str);
legend(handles.ax1,'Cell 1','Cell 2');
handles.currentcompare = 1;
[r1,c1] = find(squeeze(A(handles.r(1),:,:)));
[r2,c2] = find(squeeze(A(handles.co(1),:,:)));
axes(handles.ax2);
imshow(handles.rawdata1.results.Cn,[0 1]);
hold(handles.ax2,'on')

img = handles.rawdata1.results.Cn;
ind = [handles.r(handles.currentcompare) handles.co(handles.currentcompare)];
thr =  str2num(get(handles.thr,'String'));
with_label = false;
Coor = get_contours(handles,thr,ind,1);
handles.Coor = Coor;
d1 = handles.rawdata1.results.options.d1;
d2 = handles.rawdata1.results.options.d2;

plot_contours(handles.rawdata1.results.A(:, ind),d1,d2, img, thr,with_label, [], handles.Coor, 2);
colormap gray;
hold(handles.ax2,'off')
% x_pixel = handles.rawdata1.results.options.d1;
% y_pixel = handles.rawdata1.results.options.d2;
% xlim(handles.ax2,[0 x_pixel])
% ylim(handles.ax2,[0 y_pixel])
title(handles.ax2,'Spatial locations of the cells');



set(handles.slidertoaddtodellist,'Min',1,'Max',size(handles.updatedresults,1),'SliderStep',[1/(size(handles.updatedresults,1)-1),1],'Value',1)
str2 = sprintf('%d pairs of neurons satisfy the thresholds',length(handles.r));
set(handles.n_pair,'String',str2);

plot(handles.distcorr,pw,cc,'.r','LineStyle','none','MarkerSize',10);
hold(handles.distcorr,'on')
plot(handles.distcorr,handles.pwdist(handles.r(1),handles.co(1)),handles.crosscoef(handles.r(1),handles.co(1)),'.k','LineStyle','none','MarkerSize',15)
hold(handles.distcorr,'off')
axes(handles.histsnr)
histogram(handles.SNR);
guidata(hObject,handles)


% --- Executes on button press in addlazy.
function addlazy_Callback(hObject, eventdata, handles)
% hObject    handle to addlazy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)\
tt = handles.currentcompare;
cell_1 = handles.r(tt);
cell_2 = handles.co(tt);
m = handles.mergecells;
m = [m; {mat2str([cell_1 cell_2])} ];
set(handles.mergelist,'String',m,'Value',length(m))
handles.mergecells=m;

% add go to next pair to lazy button
C = handles.rawdata1.results.C_raw;

if handles.currentcompare~=length(handles.r)
    k = handles.currentcompare+1;
    handles.currentcompare = k;
    plot(handles.ax1,C(handles.r(k),:),'-r');
    hold(handles.ax1,'on')
    plot(handles.ax1,C(handles.co(k),:),'-b');
    hold(handles.ax1,'off')
    str = sprintf('Cell 1 - %d \n Cell 2 - %d \n Dist b/w them - %0.2d Crosscorr - %0.2d', handles.r(k),handles.co(k),handles.pwdist(handles.r(k),handles.co(k)),handles.crosscoef(handles.r(k),handles.co(k)));
    title(handles.ax1,str);
    legend(handles.ax1,'Cell 1','Cell 2');
    cc = tril(handles.crosscoef);
    pw = reshape(handles.pwdist,1,size(handles.pwdist,1)*size(handles.pwdist,2));
    cc = reshape(cc,1,size(cc,1)*size(cc,2));
    ind0 = find(pw==0);
    pw(ind0)=[];
    cc(ind0)=[];
    axes(handles.ax2);
    imshow(handles.rawdata1.results.Cn,[0 1]);
    
    img = handles.rawdata1.results.Cn;
    ind = [handles.r(handles.currentcompare) handles.co(handles.currentcompare)];
    thr =  str2num(get(handles.thr,'String'));
    with_label = false;
    Coor = get_contours(handles,thr,ind,1);
    handles.Coor = Coor;
    d1 = handles.rawdata1.results.options.d1;
    d2 = handles.rawdata1.results.options.d2;
    plot_contours(handles.updatedA(:, ind),d1,d2, img, thr,with_label, [], handles.Coor, 2);
    colormap gray;
    hold(handles.ax2,'off')
    % x_pixel = handles.rawdata1.results.options.d1;
    % y_pixel = handles.rawdata1.results.options.d2;
    % xlim(handles.ax2,[0 x_pixel])
    % ylim(handles.ax2,[0 y_pixel])
    title(handles.ax2,'Spatial locations of the cells')
    plot(handles.distcorr,pw,cc,'.r','LineStyle','none','MarkerSize',10);
    hold(handles.distcorr,'on')
    plot(handles.distcorr,handles.pwdist(handles.r(handles.currentcompare),handles.co(handles.currentcompare)),handles.crosscoef(handles.r(handles.currentcompare),handles.co(handles.currentcompare)),'.k','LineStyle','none','MarkerSize',15)
    hold(handles.distcorr,'off')
end

guidata(hObject,handles)


% --- Executes on button press in checkmerge.
function checkmerge_Callback(hObject, eventdata, handles)
% hObject    handle to checkmerge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp = [];
for i = 1:length(handles.mergecells)
    tmp = [tmp str2num(cell2mat(handles.mergecells(i)))];
end
set(handles.checkmergetxt,'String','Everything is awesome!! ');
tbl = tabulate(tmp);
if ~isempty(find(tbl(:,2)>1))
    str = sprintf('Uh oh. You have cells that are repeated in the list');
    set(handles.checkmergetxt,'String',[]);
    set(handles.checkmergetxt,'String',str);
end

guidata(hObject,handles)




function checkmergetxt_Callback(hObject, eventdata, handles)
% hObject    handle to checkmergetxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of checkmergetxt as text
%        str2double(get(hObject,'String')) returns contents of checkmergetxt as a double


% --- Executes during object creation, after setting all properties.
function checkmergetxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to checkmergetxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in multimerge.
function multimerge_Callback(hObject, eventdata, handles)
% hObject    handle to multimerge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp = [];tmp1=[];
for i = 1:length(handles.mergecells)
    tmp = [tmp str2num(cell2mat(handles.mergecells(i)))];
end
tmp1 = cell2struct(handles.mergecells,'ids',length(handles.mergecells));
tbl = tabulate(tmp);
ind_multi = tbl(find(tbl(:,2)>1),1);
keep_ind=[];
for i = 1:length(ind_multi)
    ind_m = [];
    for y = 1:size(tmp1,1)
        if ~isempty(find(str2num(tmp1(y).ids)==ind_multi(i)))
            ind_m = [ind_m y];
        end
    end
    keep_ind = [keep_ind ind_m];
    tmp_multi = [];
    for j = 1:length(ind_m)
        tmp_multi = [tmp_multi str2num(tmp1(ind_m(j)).ids)];
    end
    handles.multimergeids(i).ids= unique(tmp_multi);
end
m=handles.mergecells;
m(keep_ind)=[];
if ~isempty(keep_ind)
for k = 1:length(handles.multimergeids)
    m = [m; {mat2str(handles.multimergeids(k).ids)}];
end
end
handles.mergecells = m;
   set(handles.mergelist,'String',m,'Value',1);
guidata(hObject,handles)



% 

% --- Executes on button press in mergethemall.
function mergethemall_Callback(hObject, eventdata, handles)
% hObject    handle to mergethemall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
k1 = get(handles.id,'String');
tmp1=[];
tmp2=[];
for i = 1:length(k1)
    tmp3 = str2num(cell2mat(k1(i)));
    tmp4=[];
for u = 1:length(tmp3)
    tmp4 = [tmp4;handles.rawdata1.results.C_raw(tmp3(u),:)];
end
    merged_signal(i,:) = mean(tmp4);
    merged_contour(:,i) = handles.rawdata1.results.A(:,tmp3(1));
    tmp1 = [tmp1 tmp3(1)];
    tmp2 = [tmp2 tmp3(2:end)];
end
handles.updatedresults(tmp1,:) = merged_signal;
handles.updatedresults(tmp2,:)=[];
handles.updatedA(:,tmp1) = merged_contour;
handles.updatedA(:,tmp2)=[];
set(handles.mergelist,'String',[]);
set(handles.merge_alert,'String','Done merging! Restart the gui if you need to redo merging')
handles.pwdist=[];
handles.crosscoef=[];
handles.r=[];
handles.co=[];
guidata(hObject,handles);




% --- Executes on button press in delete4mmultimerge.
function delete4mmultimerge_Callback(hObject, eventdata, handles)
% hObject    handle to delete4mmultimerge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
k1 =(get(handles.id,'Value'));
k2 =(get(handles.id,'String'));
k2(k1,:)=[];
handles.id(k1,:)=[];
set(handles.id,'String',k2,'Value',size(k2,1));
guidata(hObject,handles)


% --- Executes on selection change in id.
function id_Callback(hObject, eventdata, handles)
% hObject    handle to id (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns id contents as cell array
%        contents{get(hObject,'Value')} returns selected item from id
k= get(handles.id,'Value')
C=  handles.rawdata1.results.C_raw;
tmp = handles.multimergeids(k).ids;
for i = 1:length(tmp)
    plot(handles.ax1,C(tmp(i),:));
    hold(handles.ax1,'on')
end
imshow(handles.rawdata1.results.Cn,[0 1]);
hold(handles.ax2,'on')
img = handles.rawdata1.results.Cn;
ind = tmp;
thr =  str2num(get(handles.thr,'String'));
with_label = false;
Coor = get_contours(handles,thr,ind,1);
handles.Coor = Coor;
d1 = handles.rawdata1.results.options.d1;
d2 = handles.rawdata1.results.options.d2;
plot_contours(handles.updatedA(:, tmp),d1,d2, img, thr,with_label, [], handles.Coor, 2);
colormap gray;

hold(handles.ax1,'off')
hold(handles.ax2,'off')

set(handles.id,'String',mat2str(tmp),'Value',1)
guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function id_CreateFcn(hObject, eventdata, handles)
% hObject    handle to id (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu6.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu6


% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton21.
function pushbutton21_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton22.
function pushbutton22_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton23.
function pushbutton23_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton24.
function pushbutton24_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton25.
function pushbutton25_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in listbox5.
function listbox5_Callback(hObject, eventdata, handles)
% hObject    handle to listbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox5


% --- Executes during object creation, after setting all properties.
function listbox5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
