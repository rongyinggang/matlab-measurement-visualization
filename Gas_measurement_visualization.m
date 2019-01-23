
function varargout = Gas_measurement_visualization(varargin)
% GAS_MEASUREMENT_VISUALIZATION MATLAB code for Gas_measurement_visualization.fig
%      GAS_MEASUREMENT_VISUALIZATION, by itself, creates a new GAS_MEASUREMENT_VISUALIZATION or raises the existing
%      singleton*.
%
%      H = GAS_MEASUREMENT_VISUALIZATION returns the handle to a new GAS_MEASUREMENT_VISUALIZATION or the handle to
%      the existing singleton*.
%
%      GAS_MEASUREMENT_VISUALIZATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GAS_MEASUREMENT_VISUALIZATION.M with the given input arguments.
%
%      GAS_MEASUREMENT_VISUALIZATION('Property','Value',...) creates a new GAS_MEASUREMENT_VISUALIZATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Gas_measurement_visualization_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Gas_measurement_visualization_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Gas_measurement_visualization

% Last Modified by GUIDE v2.5 10-Jan-2019 17:46:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Gas_measurement_visualization_OpeningFcn, ...
                   'gui_OutputFcn',  @Gas_measurement_visualization_OutputFcn, ...
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


% --- Executes just before Gas_measurement_visualization is made visible.
function Gas_measurement_visualization_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Gas_measurement_visualization (see VARARGIN)

% Choose default command line output for Gas_measurement_visualization
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Gas_measurement_visualization wait for user response (see UIRESUME)
% uiwait(handles.figure1);
 textDebug = sprintf('*** \t Push  \t "Import Button"  \t  to select the folder that all the files in the folder will be loaded automatically');    % clear debug message
 set(handles.debug, 'String', textDebug); 
 set(handles.layer_slider,'Enable','off')
 set(handles.interpolation,'Value',1)
 set(handles.data_type,'Value',1)
 set(handles.plane,'Value',1)
 
 title(handles.axes1, 'Visualization of the measurement');
 axes(handles.axes1);
 cla reset;
 
 axes(handles.axes2);
 cla reset;
 title(handles.axes2, 'Visualization of the measurement in 3D');

 view(3)
 
set(handles. figure1,'toolbar','figure');
set(handles. figure1,'menubar','figure');
% format long
global flag
flag = 0;


% --- Outputs from this function are returned to the command line.
function varargout = Gas_measurement_visualization_OutputFcn(~, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;






% --- Executes on button press in import.
function import_Callback(~, ~, handles)
% hObject    handle to import (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global flag


path = uigetdir; %gets directory
 if path==0
  % user pressed cancel
  return
 end
addpath(path);
% textDebug = sprintf('Your folder is %s', myDir);
% set(handles.debug, 'String', textDebug);

myFiles = dir(fullfile(path,'*.lvm')); %gets all csv files in struct

fileNames = {myFiles.name};
myFile = [];
% format long
global data
for k = 1:numel(fileNames) %load all the file with lvm import
    data{k} = lvm_import(fileNames{k});
    textDebug = sprintf('%s is imported',fileNames{k});
    fclose all;
    set(handles.debug, 'String', textDebug);    
end
assignin('base', 'Data', data);
assignin('base', 'FileNames', fileNames);
textDebug = sprintf('File loaded,\t  please set the parameter and press the \t "Calculation" \t button');
set(handles.debug, 'String', textDebug); 

flag = 1;

% --- Executes on button press in load_par.
function load_par_Callback(~, ~, handles)
% hObject    handle to load_par (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file,path] = uigetfile('*.txt;'); %gets directory
 if file == 0
  % user pressed cancel
  return
 end
addpath(path);

fid = fopen(file,'r');
Settings = fread(fid,'*char');
Settings = Settings';
fclose(fid);

parameter = load_parameter(Settings);
set(handles.x_abstand ,'String', parameter.schrittweite(1,1));
set(handles.y_abstand ,'String', parameter.schrittweite(1,2));
set(handles.z_abstand ,'String', parameter.schrittweite(1,3));
set(handles.x_anzahl ,'String', num2str(str2double(parameter.schritt(1,1))+1));
set(handles.y_anzahl ,'String', num2str(str2double(parameter.schritt(1,2))+1));
set(handles.z_anzahl ,'String', num2str(str2double(parameter.schritt(1,3))+1));
set(handles.Temperatur,'String',parameter.temperatur);

set(handles.a5,'String',parameter.koeffizienten(1,1));
set(handles.a4,'String',parameter.koeffizienten(1,2));
set(handles.a3,'String',parameter.koeffizienten(1,3));
set(handles.a2,'String',parameter.koeffizienten(1,4));
set(handles.a1,'String',parameter.koeffizienten(1,5));
set(handles.a0,'String',parameter.koeffizienten(1,6));
set(handles.start_x,'String',parameter.start_pos(1,1));
set(handles.start_y,'String',parameter.start_pos(1,2));
set(handles.start_z,'String',parameter.start_pos(1,3));
set(handles.sample_rate,'String',parameter.sample_rate(1,1));


textDebug = sprintf('Parameters loaded, press the \t "Calculation" \t button');
set(handles.debug, 'String', textDebug);


% --- Executes on button press in calc.
function calc_Callback(hObject, ~, handles)
% hObject    handle to calc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global flag
global data
if isempty(data)
    msgbox('Please import the files first.');
    return
    
else
    format long
    
    
    % read inputs and workspace
    fileNames = evalin('base','FileNames');
    data = evalin('base','Data');
    
    
    
    T = str2double(get(handles.Temperatur,'String'));
    a0 = str2double(get(handles.a0,'String'));
    a1 = str2double(get(handles.a1,'String'));
    a2 = str2double(get(handles.a2,'String'));
    a3 = str2double(get(handles.a3,'String'));
    a4 = str2double(get(handles.a4,'String'));
    a5 = str2double(get(handles.a5,'String'));
    x_anzahl = str2double(get(handles.x_anzahl ,'String'));
    y_anzahl  = str2double(get(handles.y_anzahl ,'String'));
    sample_rate  = str2double(get(handles.sample_rate, 'String'));
    
%     z_anzahl  = str2double(get(handles.z_anzahl ,'String'));
% 
%     x_ab = str2double(get(handles.x_abstand ,'String'));
%     y_ab = str2double(get(handles.y_abstand ,'String'));
%     z_ab = str2double(get(handles.z_abstand ,'String'));


    vel_rms = zeros(numel(fileNames),4); % define the size of the matrix

        for i = 1:numel(fileNames)
            
            len = size(data{1,i}.Segment1.data);
            
            % Temperatur value only exist in every 1000 sequence
            for j = 1:sample_rate:len(1,1)
            % fill the Temperatur value
            data{1,i}.Segment1.data(j:(j+sample_rate-1),2) = data{1,i}.Segment1.data(j,2);

            end
            
            U_kor = data{1,i}.Segment1.data(:,1).*sqrt((250-data{1,i}.Segment1.data(:,2))/(250-T));
            vel = a5 * U_kor.^5 + a4 * U_kor.^4 + a3 * U_kor.^3 + a2 * U_kor.^2 + a1 * U_kor + a0;
            vel_rms(i,1) = (rms(vel));
            vel_rms(i,2) = (mean(vel));
            vel_rms(i,3) = (std(vel));
            vel_rms(i,4) =  vel_rms(i,3)/vel_rms(i,2);
        end

    assignin('base', 'Vel_rms', vel_rms);


    n = numel(fileNames)/(x_anzahl*y_anzahl);

    h_mat = zeros(x_anzahl,y_anzahl,n);  % define the size of the matrix
    h_mat_mod = zeros(x_anzahl,y_anzahl,n);  % define the size of the matrix
    h_mat_rot = zeros(x_anzahl,y_anzahl,n);  % define the size of the matrix
    
    h_mat_m = zeros(x_anzahl,y_anzahl,n);  % define the size of the matrix
    h_mat_mod_m = zeros(x_anzahl,y_anzahl,n);  % define the size of the matrix
    h_mat_rot_m = zeros(x_anzahl,y_anzahl,n);  % define the size of the matrix
    
    h_mat_std = zeros(x_anzahl,y_anzahl,n);  % define the size of the matrix
    h_mat_mod_std = zeros(x_anzahl,y_anzahl,n);  % define the size of the matrix
    h_mat_rot_std = zeros(x_anzahl,y_anzahl,n);  % define the size of the matrix
    
    h_mat_tb = zeros(x_anzahl,y_anzahl,n);  % define the size of the matrix
    h_mat_mod_tb = zeros(x_anzahl,y_anzahl,n);  % define the size of the matrix
    h_mat_rot_tb = zeros(x_anzahl,y_anzahl,n);  % define the size of the matrix
    
    
    %% rms
    
    for j=0:n-1

    h_mat(:,:,j+1) = reshape(vel_rms((j*x_anzahl*y_anzahl)+1:x_anzahl*y_anzahl*(j+1),1),[x_anzahl,y_anzahl]);

    end

    % modification of the matrix to fit the measurement route

    for k=0:n-1

    h_mat_mod(:,1:2:end,k+1) = h_mat(:,1:2:end,k+1);

    h_mat_mod(:,2:2:end,k+1) = flipud(h_mat(:,2:2:end,k+1));

    end

    for l=0:n-1

    h_mat_rot(:,:,l+1) = rot90(h_mat_mod(:,:,l+1),1);

    end
    
    %% mean
    for j=0:n-1

    h_mat_m(:,:,j+1) = reshape(vel_rms((j*x_anzahl*y_anzahl)+1:x_anzahl*y_anzahl*(j+1),2),[x_anzahl,y_anzahl]);

    end


    % modification of the matrix to fit the measurement route

    for k=0:n-1

    h_mat_mod_m(:,1:2:end,k+1) = h_mat_m(:,1:2:end,k+1);

    h_mat_mod_m(:,2:2:end,k+1) = flipud(h_mat_m(:,2:2:end,k+1));

    end

    for l=0:n-1

    h_mat_rot_m(:,:,l+1) = rot90(h_mat_mod_m(:,:,l+1),1);

    end
    
    %% STD
    for j=0:n-1

    h_mat_std(:,:,j+1) = reshape(vel_rms((j*x_anzahl*y_anzahl)+1:x_anzahl*y_anzahl*(j+1),3),[x_anzahl,y_anzahl]);

    end

    % modification of the matrix to fit the measurement route

    for k=0:n-1

    h_mat_mod_std(:,1:2:end,k+1) = h_mat_std(:,1:2:end,k+1);

    h_mat_mod_std(:,2:2:end,k+1) = flipud(h_mat_std(:,2:2:end,k+1));

    end

    for l=0:n-1

    h_mat_rot_std(:,:,l+1) = rot90(h_mat_mod_std(:,:,l+1),1);

    end
    
    %% Turbulence 
    
    for j=0:n-1

    h_mat_tb(:,:,j+1) = reshape(vel_rms((j*x_anzahl*y_anzahl)+1:x_anzahl*y_anzahl*(j+1),4),[x_anzahl,y_anzahl]);

    end

    % modification of the matrix to fit the measurement route

    for k=0:n-1

    h_mat_mod_tb(:,1:2:end,k+1) = h_mat_tb(:,1:2:end,k+1);

    h_mat_mod_tb(:,2:2:end,k+1) = flipud(h_mat_tb(:,2:2:end,k+1));

    end

    for l=0:n-1

    h_mat_rot_tb(:,:,l+1) = rot90(h_mat_mod_tb(:,:,l+1),1);

    end
    %% pass the value to workspace 

%     assignin('base', 'H_mat', h_mat);
%     assignin('base', 'H_mat_rot', h_mat_rot);

    h_mat_xy = flipud(h_mat_rot);
    h_mat_xz = permute(h_mat_xy,[3 2 1]); % rotate the matrix fit the correct coordinate
    h_mat_yz = permute(h_mat_xy,[3 1 2]); % rotate the matrix to fit the correct coordinate
    
    h_mat_xy_m = flipud(h_mat_rot_m);
    h_mat_xz_m = permute(h_mat_xy_m,[3 2 1]); % rotate the matrix fit the correct coordinate
    h_mat_yz_m = permute(h_mat_xy_m,[3 1 2]); % rotate the matrix to fit the correct coordinate

    h_mat_xy_std = flipud(h_mat_rot_std);
    h_mat_xz_std = permute(h_mat_xy_std,[3 2 1]); % rotate the matrix fit the correct coordinate
    h_mat_yz_std = permute(h_mat_xy_std,[3 1 2]); % rotate the matrix to fit the correct coordinate

    h_mat_xy_tb = flipud(h_mat_rot_tb);
    h_mat_xz_tb = permute(h_mat_xy_tb,[3 2 1]); % rotate the matrix fit the correct coordinate
    h_mat_yz_tb = permute(h_mat_xy_tb,[3 1 2]); % rotate the matrix to fit the correct coordinate
    
    assignin('base', 'H_mat_xy', h_mat_xy);
    assignin('base', 'H_mat_xz', h_mat_xz);
    assignin('base', 'H_mat_yz', h_mat_yz);
    
    assignin('base', 'H_mat_xy_m', h_mat_xy_m);
    assignin('base', 'H_mat_xz_m', h_mat_xz_m);
    assignin('base', 'H_mat_yz_m', h_mat_yz_m);
    
    assignin('base', 'H_mat_xy_std', h_mat_xy_std);
    assignin('base', 'H_mat_xz_std', h_mat_xz_std);
    assignin('base', 'H_mat_yz_std', h_mat_yz_std);
    
    assignin('base', 'H_mat_xy_tb', h_mat_xy_tb);
    assignin('base', 'H_mat_xz_tb', h_mat_xz_tb);
    assignin('base', 'H_mat_yz_tb', h_mat_yz_tb);
    

    textDebug = sprintf('Calculation finished, now you can choose a Plane and a Layer to see the results');
    set(handles.debug, 'String', textDebug); 
    flag = 2;

end


% --- Executes on selection change in plane.
function plane_Callback(hObject, eventdata, handles)
% hObject    handle to plane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plane contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plane

global flag
if flag == 2
    
    z_start = str2double(get(handles.start_z ,'String'));
    
    x_anzahl = str2double(get(handles.x_anzahl ,'String'));
    y_anzahl  = str2double(get(handles.y_anzahl ,'String'));
    z_anzahl  = str2double(get(handles.z_anzahl ,'String'));

    x_ab = str2double(get(handles.x_abstand ,'String'));
    y_ab = str2double(get(handles.y_abstand ,'String'));
    z_ab = str2double(get(handles.z_abstand ,'String'));
    
    

    x_limit_sb = [1 x_anzahl];
    y_limit_sb = [1 y_anzahl];
    z_limit_sb = [1 z_anzahl];

    x_limit = [-(x_anzahl-1)*x_ab/2 (x_anzahl-1)*x_ab/2];
    y_limit = [-(y_anzahl-1)*y_ab/2 (y_anzahl-1)*y_ab/2];
    z_limit = [z_start z_start + (z_anzahl-1)*z_ab];
    
    x_tick = linspace(-(x_anzahl-1)*x_ab/2,(x_anzahl-1)*x_ab/2,4);
    y_tick = linspace(-(y_anzahl-1)*y_ab/2,(y_anzahl-1)*y_ab/2,4);
    z_tick = linspace(z_start, z_start + (z_anzahl-1)*z_ab,3);

    h_mat_xy = evalin('base','H_mat_xy');
    h_mat_xz = evalin('base','H_mat_xz');
    h_mat_yz = evalin('base','H_mat_yz');
    h_mat_xy_m = evalin('base','H_mat_xy_m');
    h_mat_xz_m = evalin('base','H_mat_xz_m');
    h_mat_yz_m = evalin('base','H_mat_yz_m');
    h_mat_xy_std = evalin('base','H_mat_xy_std');
    h_mat_xz_std = evalin('base','H_mat_xz_std');
    h_mat_yz_std = evalin('base','H_mat_yz_std');
    h_mat_xy_tb = evalin('base','H_mat_xy_tb');
    h_mat_xz_tb = evalin('base','H_mat_xz_tb');
    h_mat_yz_tb = evalin('base','H_mat_yz_tb');
    vel_rms = evalin('base','Vel_rms');

    interpolation =  get(handles.interpolation,'Value');
    plane = get(handles.plane,'Value');
    data_type = get(handles.data_type,'Value');
    
    x = -(x_anzahl-1)/2:(x_anzahl-1)/2;
    y = -(y_anzahl-1)/2:(y_anzahl-1)/2;
    z = 0 :(z_anzahl-1);  
    
     switch data_type
        case 1
             colorrange = [min(vel_rms(:,1)) max(vel_rms(:,1))];
             Colorbar_label = sprintf('Root Mean Square Velocity [m/s]');
             title_label =  sprintf('Root Mean Square Velocity ');
  
        case 2
             colorrange = [min(vel_rms(:,2)) max(vel_rms(:,2))];
             h_mat_xy = h_mat_xy_m;
             h_mat_xz = h_mat_xz_m;
             h_mat_yz = h_mat_yz_m;
             Colorbar_label = sprintf('Mean Velocity [m/s]');
             title_label =  sprintf('Mean Velocity ');
        case 3
             colorrange = [min(vel_rms(:,3)) max(vel_rms(:,3))];
             h_mat_xy = h_mat_xy_std;
             h_mat_xz = h_mat_xz_std;
             h_mat_yz = h_mat_yz_std;
             Colorbar_label = sprintf('Standard Deviation [m/s]');
             title_label =  sprintf('Standard Deviation ');
             
        case 4
             colorrange = [min(vel_rms(:,4)) max(vel_rms(:,4))];
             h_mat_xy = h_mat_xy_tb;
             h_mat_xz = h_mat_xz_tb;
             h_mat_yz = h_mat_yz_tb;
             Colorbar_label = sprintf('Turbulence Grad ');
             title_label =  sprintf('Turbulence Grad ');
    
    end

    switch plane

        case 1
            
            textDebug = sprintf('Please select a plane');
            set(handles.debug, 'String', textDebug); 
            set(handles.layer_slider,'Enable','off')
            set(handles.text_layer,'String', '');

        case 2
            textDebug = sprintf('');
            set(handles.debug, 'String', textDebug); 

            set(handles.layer_slider,'Enable','on')
            set(handles.layer_slider, 'min', z_limit_sb(1,1));
            set(handles.layer_slider, 'max', z_limit_sb(1,2));
            set(handles.layer_slider, 'Value', z_limit_sb(1,1));                      
            set(handles.layer_slider, 'SliderStep',[1/(z_anzahl-1),1/(z_anzahl-1)]);           
            
            %% Plot 2D figure  
            axes(handles.axes1);
            textDebug = sprintf('');
            set(handles.debug, 'String', textDebug); 
            value = get(hObject,'Value');                  
            value = round(value);                                 
            set(hObject, 'Value', value);  
            z_layer =  get(handles.layer_slider,'value');                                            
            text_slider = sprintf('[ Z = %d mm ]',z_start + (z_layer-1)*z_ab);
            set(handles.text_layer,'String', text_slider); 
            [X,Y] = meshgrid(x*x_ab,y*y_ab);
            Figure = surf(X,Y,h_mat_xy(:,:,z_layer));
            text_title = strcat(title_label,' [ Z = ',' ',num2str(z_start + (z_layer-1)*z_ab),' mm ]');
            title(text_title);
            xlabel('x direction [mm]') 
            ylabel('y direction [mm]') 
            set(handles.axes1,'Xlim',x_limit,'Ylim',y_limit,'XTick',x_tick,'YTick',y_tick);
            view(2) % view 2D surf plot

            if interpolation == 1 
                shading interp
            else
                shading flat
            end
            
            colormap(jet)
            caxis(colorrange)
            

            %% Plot 3D figure 
            axes(handles.axes2);
            [X,Y,Z] = meshgrid(x*x_ab,y*y_ab,z*z_ab + z_start);
            slice(X,Y,Z,h_mat_xy,[],[],z_start + (z_layer-1)*z_ab)
            set(handles.axes2,'Xlim',x_limit,'Ylim',y_limit,'Zlim',z_limit,'XTick',x_tick,'YTick',y_tick,'ZTick',z_tick)
            
            title(title_label);
            xlabel('x direction [mm]','Rotation',17) 
            ylabel('y direction [mm]','Rotation',-28) 
            zlabel('z direction [mm]') 
            
            colormap(jet)
            ylabel(colorbar, Colorbar_label)
            caxis(colorrange)
            
            if interpolation == 1 
                shading interp
            else
                shading flat
            end
            


        case 3

            set(handles.layer_slider,'Enable','on')
            textDebug = sprintf('');
            set(handles.debug, 'String', textDebug); 

            set(handles.layer_slider,'Enable','on')
            set(handles.layer_slider, 'min', y_limit_sb(1,1));
            set(handles.layer_slider, 'max', y_limit_sb(1,2));
            set(handles.layer_slider, 'Value', y_limit_sb(1,1));                      
            set(handles.layer_slider, 'SliderStep',[1/(y_anzahl-1),1/(y_anzahl-1)]);           
            
            %% Plot 2D figure
            axes(handles.axes1);
            % set (gcf, 'WindowButtonMotionFcn', @mouseMove);
            textDebug = sprintf('');
            set(handles.debug, 'String', textDebug); 

            value = get(hObject,'Value');                  
            value = round(value);                                 
            set(hObject, 'Value', value);  
            y_layer =  get(handles.layer_slider,'value');                                            
            text_slider = sprintf('[ Y = %d mm ]',(get(handles.layer_slider,'value')-(y_anzahl+1)/2)*y_ab);
            set(handles.text_layer,'String', text_slider); 

            [X,Z] = meshgrid(x*x_ab,z_start + z*z_ab);
            Figure = surf(X,Z,h_mat_xz(:,:,y_layer));
            
            text_title = strcat(title_label,' [ Y = ',' ',num2str((get(handles.layer_slider,'value')-(y_anzahl+1)/2)*y_ab),' mm ]');
            title(text_title)
            xlabel('x direction [mm]') 
            ylabel('z direction [mm]') 
            set(handles.axes1,'Xlim',x_limit,'Ylim',z_limit,'XTick',x_tick,'YTick',z_tick);
            view(2) % view 2D surf plot

            if interpolation == 1 
                shading interp
            else
                shading flat
            end
            
            colormap(jet)
            caxis(colorrange)
            
            %% Plot 3D figure 
            axes(handles.axes2);
            [X,Y,Z] = meshgrid(x*x_ab,y*y_ab,z_start + z*z_ab);
            slice(X,Y,Z,h_mat_xy,[],((get(handles.layer_slider,'value')-(y_anzahl+1)/2)*y_ab),[])
            set(handles.axes2,'Xlim',x_limit,'Ylim',y_limit,'Zlim',z_limit,'XTick',x_tick,'YTick',y_tick,'ZTick',z_tick)
            
            title(title_label);
            xlabel('x direction [mm]','Rotation',17) 
            ylabel('y direction [mm]','Rotation',-28) 
            zlabel('z direction [mm]') 
            
            colormap(jet)
            ylabel(colorbar, Colorbar_label)
            caxis(colorrange)
            
            if interpolation == 1 
                shading interp
            else
                shading flat
            end



        case 4

            set(handles.layer_slider,'Enable','on')
            textDebug = sprintf('');
            set(handles.debug, 'String', textDebug); 

            set(handles.layer_slider,'Enable','on')
            set(handles.layer_slider, 'min', x_limit_sb(1,1));
            set(handles.layer_slider, 'max', x_limit_sb(1,2));
            set(handles.layer_slider, 'Value', x_limit_sb(1,1));                      
            set(handles.layer_slider, 'SliderStep',[1/(x_anzahl-1),1/(x_anzahl-1)]);           
            
            %% Plot 2d figure
            axes(handles.axes1);
            % set (gcf, 'WindowButtonMotionFcn', @mouseMove);
            textDebug = sprintf('');
            set(handles.debug, 'String', textDebug); 

            value = get(hObject,'Value');                  
            value = round(value);                                 
            set(hObject, 'Value', value);  
            x_layer =  get(handles.layer_slider,'value');                                            
            text_slider = sprintf('[ X = %d mm ]',(get(handles.layer_slider,'value')-(x_anzahl+1)/2)*x_ab);
            set(handles.text_layer,'String', text_slider); 
            [Y,Z] = meshgrid(y*y_ab,z_start + z*z_ab);
            Figure = surf(Y,Z,h_mat_yz(:,:,x_layer));
            text_title = strcat(title_label,' [ X = ',' ',num2str((get(handles.layer_slider,'value')-(x_anzahl+1)/2)*x_ab),' mm ]');
            title(text_title);
            xlabel('y direction [mm]') 
            ylabel('z direction [mm]')
            set(handles.axes1,'Xlim',y_limit,'Ylim',z_limit,'XTick',y_tick,'YTick',z_tick);
            view(2) % view 2D surf plot     

            if interpolation == 1 
                shading interp
            else
                shading flat
            end
            
            colormap(jet)
            caxis(colorrange)
            
            %% Plot 3D figure 
            axes(handles.axes2);
            [X,Y,Z] = meshgrid(x*x_ab,y*y_ab,z_start + z*z_ab);
            slice(X,Y,Z,h_mat_xy,((get(handles.layer_slider,'value')-(x_anzahl+1)/2)*x_ab),[],[])
            set(handles.axes2,'Xlim',x_limit,'Ylim',y_limit,'Zlim',z_limit,'XTick',x_tick,'YTick',y_tick,'ZTick',z_tick)
            
            title(title_label);
            xlabel('x direction [mm]','Rotation',17) 
            ylabel('y direction [mm]','Rotation',-28)  
            zlabel('z direction [mm]') 
            
            colormap(jet)
            ylabel(colorbar, Colorbar_label)
            caxis(colorrange)
            
            if interpolation == 1 
                shading interp
            else
                shading flat
            end
    end
else 
    msgbox('Please import the files and calculate the results first.');
    return
end

    
    
    
    



% --- Executes on slider movement.
function layer_slider_Callback(hObject, eventdata, handles)
% hObject    handle to layer_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global flag
if flag == 2
    
    z_start =  str2double(get(handles.start_z ,'String'));
    
    x_anzahl = str2double(get(handles.x_anzahl ,'String'));
    y_anzahl  = str2double(get(handles.y_anzahl ,'String'));
    z_anzahl  = str2double(get(handles.z_anzahl ,'String'));

    x_ab = str2double(get(handles.x_abstand ,'String'));
    y_ab = str2double(get(handles.y_abstand ,'String'));
    z_ab = str2double(get(handles.z_abstand ,'String'));

    x_limit_sb = [1 x_anzahl];
    y_limit_sb = [1 y_anzahl];
    z_limit_sb = [1 z_anzahl];

    x_limit = [-(x_anzahl-1)*x_ab/2 (x_anzahl-1)*x_ab/2];
    y_limit = [-(y_anzahl-1)*y_ab/2 (y_anzahl-1)*y_ab/2];
    z_limit = [z_start z_start + (z_anzahl-1)*z_ab];
    
    x_tick = linspace(-(x_anzahl-1)*x_ab/2,(x_anzahl-1)*x_ab/2,4);
    y_tick = linspace(-(y_anzahl-1)*y_ab/2,(y_anzahl-1)*y_ab/2,4);
    z_tick = linspace(z_start, z_start + (z_anzahl-1)*z_ab,3);
    
    h_mat_xy = evalin('base','H_mat_xy');
    h_mat_xz = evalin('base','H_mat_xz');
    h_mat_yz = evalin('base','H_mat_yz');
    h_mat_xy_m = evalin('base','H_mat_xy_m');
    h_mat_xz_m = evalin('base','H_mat_xz_m');
    h_mat_yz_m = evalin('base','H_mat_yz_m');
    h_mat_xy_std = evalin('base','H_mat_xy_std');
    h_mat_xz_std = evalin('base','H_mat_xz_std');
    h_mat_yz_std = evalin('base','H_mat_yz_std');
    h_mat_xy_tb = evalin('base','H_mat_xy_tb');
    h_mat_xz_tb = evalin('base','H_mat_xz_tb');
    h_mat_yz_tb = evalin('base','H_mat_yz_tb');
    vel_rms = evalin('base','Vel_rms');

    interpolation =  get(handles.interpolation,'Value');
    plane = get(handles.plane,'Value');
    data_type = get(handles.data_type,'Value');
    
    x = -(x_anzahl-1)/2:(x_anzahl-1)/2;
    y = -(y_anzahl-1)/2:(y_anzahl-1)/2;
    z = 0 :(z_anzahl-1);  
    
     switch data_type
        case 1
             colorrange = [min(vel_rms(:,1)) max(vel_rms(:,1))];
             Colorbar_label = sprintf('Root Mean Square Velocity [m/s]');
             title_label =  sprintf('Root Mean Square Velocity ');
  
        case 2
             colorrange = [min(vel_rms(:,2)) max(vel_rms(:,2))];
             h_mat_xy = h_mat_xy_m;
             h_mat_xz = h_mat_xz_m;
             h_mat_yz = h_mat_yz_m;
             Colorbar_label = sprintf('Mean Velocity [m/s]');
             title_label =  sprintf('Mean Velocity ');
        case 3
             colorrange = [min(vel_rms(:,3)) max(vel_rms(:,3))];
             h_mat_xy = h_mat_xy_std;
             h_mat_xz = h_mat_xz_std;
             h_mat_yz = h_mat_yz_std;
             Colorbar_label = sprintf('Standard Deviation [m/s]');
             title_label =  sprintf('Standard Deviation ');
        case 4
             colorrange = [min(vel_rms(:,4)) max(vel_rms(:,4))];
             h_mat_xy = h_mat_xy_tb;
             h_mat_xz = h_mat_xz_tb;
             h_mat_yz = h_mat_yz_tb;
             Colorbar_label = sprintf('Turbulence Grad ');
             title_label =  sprintf('Turbulence Grad ');
    
    end

    switch plane

        case 1
            
            textDebug = sprintf('Please select a plane');
            set(handles.debug, 'String', textDebug); 
            set(handles.layer_slider,'Enable','off')
            set(handles.text_layer,'String', '');

        case 2
            textDebug = sprintf('');
            set(handles.debug, 'String', textDebug); 

            set(handles.layer_slider,'Enable','on')
            set(handles.layer_slider, 'min', z_limit_sb(1,1));
            set(handles.layer_slider, 'max', z_limit_sb(1,2));
                                 
            set(handles.layer_slider, 'SliderStep',[1/(z_anzahl-1),1/(z_anzahl-1)]);           
            
            %% Plot 2D figure  
            axes(handles.axes1);
            textDebug = sprintf('');
            set(handles.debug, 'String', textDebug); 
            value = get(hObject,'Value');                  
            value = round(value);                                 
            set(hObject, 'Value', value);  
            z_layer =  get(handles.layer_slider,'value');                                            
            text_slider = sprintf('[ Z = %d mm ]',z_start + (z_layer-1)*z_ab);
            % text_slider = sprintf('[ Z = %d mm ]',(get(handles.layer_slider,'value')*z_ab + z_start));
            set(handles.text_layer,'String', text_slider); 
            [X,Y] = meshgrid(x*x_ab,y*y_ab);
            Figure = surf(X,Y,h_mat_xy(:,:,z_layer));
            text_title = strcat(title_label,' [ Z = ',' ',num2str(z_start + (z_layer-1)*z_ab),' mm ]');
  
            title(text_title);
            xlabel('x direction [mm]') 
            ylabel('y direction [mm]') 
            set(handles.axes1,'Xlim',x_limit,'Ylim',y_limit,'XTick',x_tick,'YTick',y_tick);
            view(2) % view 2D surf plot

            if interpolation == 1 
                shading interp
            else
                shading flat
            end
            
            colormap(jet)
            caxis(colorrange)
            

            %% Plot 3D figure 
            axes(handles.axes2);
            [X,Y,Z] = meshgrid(x*x_ab,y*y_ab,z_start + z*z_ab);
            slice(X,Y,Z,h_mat_xy,[],[],z_start + (z_layer-1)*z_ab)
            set(handles.axes2,'Xlim',x_limit,'Ylim',y_limit,'Zlim',z_limit,'XTick',x_tick,'YTick',y_tick,'ZTick',z_tick)
            
            title(title_label);
            xlabel('x direction [mm]','Rotation',17) 
            ylabel('y direction [mm]','Rotation',-28)  
            zlabel('z direction [mm]') 
            
            colormap(jet)
            ylabel(colorbar, Colorbar_label)
            caxis(colorrange)
            
            if interpolation == 1 
                shading interp
            else
                shading flat
            end
            


        case 3

            set(handles.layer_slider,'Enable','on')
            textDebug = sprintf('');
            set(handles.debug, 'String', textDebug); 

            set(handles.layer_slider,'Enable','on')
            set(handles.layer_slider, 'min', y_limit_sb(1,1));
            set(handles.layer_slider, 'max', y_limit_sb(1,2));
                                 
            set(handles.layer_slider, 'SliderStep',[1/(y_anzahl-1),1/(y_anzahl-1)]);           
            
            %% Plot 2D figure
            axes(handles.axes1);
            % set (gcf, 'WindowButtonMotionFcn', @mouseMove);
            textDebug = sprintf('');
            set(handles.debug, 'String', textDebug); 

            value = get(hObject,'Value');                  
            value = round(value);                                 
            set(hObject, 'Value', value);  
            y_layer =  get(handles.layer_slider,'value');                                            
            text_slider = sprintf('[ Y = %d mm ]',(get(handles.layer_slider,'value')-(y_anzahl+1)/2)*y_ab);
            set(handles.text_layer,'String', text_slider); 

            [X,Z] = meshgrid(x*x_ab,z_start + z*z_ab);
            Figure = surf(X,Z,h_mat_xz(:,:,y_layer));
            
            text_title = strcat(title_label,' [ Y = ',' ',num2str((get(handles.layer_slider,'value')-(y_anzahl+1)/2)*y_ab),' mm ]');
            title(text_title)
            xlabel('x direction [mm]') 
            ylabel('z direction [mm]') 
            set(handles.axes1,'Xlim',x_limit,'Ylim',z_limit,'XTick',x_tick,'YTick',z_tick);
            view(2) % view 2D surf plot

            if interpolation == 1 
                shading interp
            else
                shading flat
            end
            
            colormap(jet)
            caxis(colorrange)
            
            %% Plot 3D figure 
            axes(handles.axes2);
            [X,Y,Z] = meshgrid(x*x_ab,y*y_ab,z_start + z*z_ab);
            slice(X,Y,Z,h_mat_xy,[],((get(handles.layer_slider,'value')-(y_anzahl+1)/2)*y_ab),[])
            set(handles.axes2,'Xlim',x_limit,'Ylim',y_limit,'Zlim',z_limit,'XTick',x_tick,'YTick',y_tick,'ZTick',z_tick)
            
            title(title_label);
            xlabel('x direction [mm]','Rotation',17) 
            ylabel('y direction [mm]','Rotation',-28) 
            zlabel('z direction [mm]') 
            
            colormap(jet)
            ylabel(colorbar, Colorbar_label)
            caxis(colorrange)
            
            if interpolation == 1 
                shading interp
            else
                shading flat
            end



        case 4

            set(handles.layer_slider,'Enable','on')
            textDebug = sprintf('');
            set(handles.debug, 'String', textDebug); 

            set(handles.layer_slider,'Enable','on')
            set(handles.layer_slider, 'min', x_limit_sb(1,1));
            set(handles.layer_slider, 'max', x_limit_sb(1,2));
                                 
            set(handles.layer_slider, 'SliderStep',[1/(x_anzahl-1),1/(x_anzahl-1)]);           
            
            %% Plot 2d figure
            axes(handles.axes1);
            % set (gcf, 'WindowButtonMotionFcn', @mouseMove);
            textDebug = sprintf('');
            set(handles.debug, 'String', textDebug); 

            value = get(hObject,'Value');                  
            value = round(value);                                 
            set(hObject, 'Value', value);  
            x_layer =  get(handles.layer_slider,'value');                                            
            text_slider = sprintf('[ X = %d mm ]',(get(handles.layer_slider,'value')-(x_anzahl+1)/2)*x_ab);
            set(handles.text_layer,'String', text_slider); 
            [Y,Z] = meshgrid(y*y_ab,z_start + z*z_ab);
            Figure = surf(Y,Z,h_mat_yz(:,:,x_layer));
            text_title = strcat(title_label,' [ X = ',' ',num2str((get(handles.layer_slider,'value')-(x_anzahl+1)/2)*x_ab),' mm ]');
            title(text_title);
            xlabel('y direction [mm]') 
            ylabel('z direction [mm]')
            set(handles.axes1,'Xlim',y_limit,'Ylim',z_limit,'XTick',y_tick,'YTick',z_tick);
            view(2) % view 2D surf plot     

            if interpolation == 1 
                shading interp
            else
                shading flat
            end
            
            colormap(jet)
            caxis(colorrange)
            
            %% Plot 3D figure 
            axes(handles.axes2);
            [X,Y,Z] = meshgrid(x*x_ab,y*y_ab,z_start + z*z_ab);
            slice(X,Y,Z,h_mat_xy,((get(handles.layer_slider,'value')-(x_anzahl+1)/2)*x_ab),[],[])
            set(handles.axes2,'Xlim',x_limit,'Ylim',y_limit,'Zlim',z_limit,'XTick',x_tick,'YTick',y_tick,'ZTick',z_tick)
            
            title(title_label);
            xlabel('x direction [mm]','Rotation',17) 
            ylabel('y direction [mm]','Rotation',-28) 
            zlabel('z direction [mm]') 
            
            colormap(jet)
            ylabel(colorbar, Colorbar_label)
            caxis(colorrange)
            
            if interpolation == 1 
                shading interp
            else
                shading flat
            end
    end
else 
    msgbox('Please import the files and calculate the results first.');
    return
end








% --- Executes on selection change in interpolation.
function interpolation_Callback(hObject, eventdata, handles)
% hObject    handle to interpolation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns interpolation contents as cell array
%        contents{get(hObject,'Value')} returns selected item from interpolation
global flag
if flag == 2
    z_start =  str2double(get(handles.start_z ,'String'));
    x_anzahl = str2double(get(handles.x_anzahl ,'String'));
    y_anzahl  = str2double(get(handles.y_anzahl ,'String'));
    z_anzahl  = str2double(get(handles.z_anzahl ,'String'));

    x_ab = str2double(get(handles.x_abstand ,'String'));
    y_ab = str2double(get(handles.y_abstand ,'String'));
    z_ab = str2double(get(handles.z_abstand ,'String'));

    x_limit_sb = [1 x_anzahl];
    y_limit_sb = [1 y_anzahl];
    z_limit_sb = [1 z_anzahl];

    x_limit = [-(x_anzahl-1)*x_ab/2 (x_anzahl-1)*x_ab/2];
    y_limit = [-(y_anzahl-1)*y_ab/2 (y_anzahl-1)*y_ab/2];
    z_limit = [z_start z_start + (z_anzahl-1)*z_ab];
    
    x_tick = linspace(-(x_anzahl-1)*x_ab/2,(x_anzahl-1)*x_ab/2,4);
    y_tick = linspace(-(y_anzahl-1)*y_ab/2,(y_anzahl-1)*y_ab/2,4);
    z_tick = linspace(z_start, z_start + (z_anzahl-1)*z_ab,3);

    h_mat_xy = evalin('base','H_mat_xy');
    h_mat_xz = evalin('base','H_mat_xz');
    h_mat_yz = evalin('base','H_mat_yz');
    h_mat_xy_m = evalin('base','H_mat_xy_m');
    h_mat_xz_m = evalin('base','H_mat_xz_m');
    h_mat_yz_m = evalin('base','H_mat_yz_m');
    h_mat_xy_std = evalin('base','H_mat_xy_std');
    h_mat_xz_std = evalin('base','H_mat_xz_std');
    h_mat_yz_std = evalin('base','H_mat_yz_std');
    h_mat_xy_tb = evalin('base','H_mat_xy_tb');
    h_mat_xz_tb = evalin('base','H_mat_xz_tb');
    h_mat_yz_tb = evalin('base','H_mat_yz_tb');
    vel_rms = evalin('base','Vel_rms');

    interpolation =  get(handles.interpolation,'Value');
    plane = get(handles.plane,'Value');
    data_type = get(handles.data_type,'Value');
    
    x = -(x_anzahl-1)/2:(x_anzahl-1)/2;
    y = -(y_anzahl-1)/2:(y_anzahl-1)/2;
    z = 0 :(z_anzahl-1);  
    
     switch data_type
        case 1
             colorrange = [min(vel_rms(:,1)) max(vel_rms(:,1))];
             Colorbar_label = sprintf('Root Mean Square Velocity [m/s]');
             title_label =  sprintf('Root Mean Square Velocity ');
  
        case 2
             colorrange = [min(vel_rms(:,2)) max(vel_rms(:,2))];
             h_mat_xy = h_mat_xy_m;
             h_mat_xz = h_mat_xz_m;
             h_mat_yz = h_mat_yz_m;
             Colorbar_label = sprintf('Mean Velocity [m/s]');
             title_label =  sprintf('Mean Velocity ');
        case 3
             colorrange = [min(vel_rms(:,3)) max(vel_rms(:,3))];
             h_mat_xy = h_mat_xy_std;
             h_mat_xz = h_mat_xz_std;
             h_mat_yz = h_mat_yz_std;
             Colorbar_label = sprintf('Standard Deviation [m/s]');
             title_label =  sprintf('Standard Deviation ');
             
        case 4
             colorrange = [min(vel_rms(:,4)) max(vel_rms(:,4))];
             h_mat_xy = h_mat_xy_tb;
             h_mat_xz = h_mat_xz_tb;
             h_mat_yz = h_mat_yz_tb;
             Colorbar_label = sprintf('Turbulence Grad ');
             title_label =  sprintf('Turbulence Grad ');
    
    end

    switch plane

        case 1
            
            textDebug = sprintf('Please select a plane');
            set(handles.debug, 'String', textDebug); 
            set(handles.layer_slider,'Enable','off')
            set(handles.text_layer,'String', '');

        case 2
            textDebug = sprintf('');
            set(handles.debug, 'String', textDebug); 

            set(handles.layer_slider,'Enable','on')
            set(handles.layer_slider, 'min', z_limit_sb(1,1));
            set(handles.layer_slider, 'max', z_limit_sb(1,2));
                                
            set(handles.layer_slider, 'SliderStep',[1/(z_anzahl-1),1/(z_anzahl-1)]);           
            
            %% Plot 2D figure  
            axes(handles.axes1);
            textDebug = sprintf('');
            set(handles.debug, 'String', textDebug); 
            value = get(hObject,'Value');                  
            value = round(value);                                 
            set(hObject, 'Value', value);  
            z_layer =  get(handles.layer_slider,'value');                                            
            text_slider = sprintf('[ Z = %d mm ]',z_start + (z_layer-1)*z_ab);
            set(handles.text_layer,'String', text_slider); 
            [X,Y] = meshgrid(x*x_ab,y*y_ab);
            Figure = surf(X,Y,h_mat_xy(:,:,z_layer));
            text_title = strcat(title_label,' [ Z = ',' ',num2str(z_start + (z_layer-1)*z_ab),' mm ]');
            title(text_title);
            xlabel('x direction [mm]') 
            ylabel('y direction [mm]') 
            set(handles.axes1,'Xlim',x_limit,'Ylim',y_limit,'XTick',x_tick,'YTick',y_tick);
            view(2) % view 2D surf plot

            if interpolation == 1 
                shading interp
            else
                shading flat
            end
            
            colormap(jet)
            caxis(colorrange)
            

            %% Plot 3D figure 
            axes(handles.axes2);
            [X,Y,Z] = meshgrid(x*x_ab,y*y_ab,z_start + z*z_ab);
            slice(X,Y,Z,h_mat_xy,[],[],z_start + (z_layer-1)*z_ab)
            set(handles.axes2,'Xlim',x_limit,'Ylim',y_limit,'Zlim',z_limit,'XTick',x_tick,'YTick',y_tick,'ZTick',z_tick)
            
            title(title_label);
            xlabel('x direction [mm]','Rotation',17) 
            ylabel('y direction [mm]','Rotation',-28) 
            zlabel('z direction [mm]') 
            
            colormap(jet)
            ylabel(colorbar, Colorbar_label)
            caxis(colorrange)
            
            if interpolation == 1 
                shading interp
            else
                shading flat
            end
            


        case 3

            set(handles.layer_slider,'Enable','on')
            textDebug = sprintf('');
            set(handles.debug, 'String', textDebug); 

            set(handles.layer_slider,'Enable','on')
            set(handles.layer_slider, 'min', y_limit_sb(1,1));
            set(handles.layer_slider, 'max', y_limit_sb(1,2));
                                
            set(handles.layer_slider, 'SliderStep',[1/(y_anzahl-1),1/(y_anzahl-1)]);           
            
            %% Plot 2D figure
            axes(handles.axes1);
            % set (gcf, 'WindowButtonMotionFcn', @mouseMove);
            textDebug = sprintf('');
            set(handles.debug, 'String', textDebug); 

            value = get(hObject,'Value');                  
            value = round(value);                                 
            set(hObject, 'Value', value);  
            y_layer =  get(handles.layer_slider,'value');                                            
            text_slider = sprintf('[ Y = %d mm ]',(get(handles.layer_slider,'value')-(y_anzahl+1)/2)*y_ab);
            set(handles.text_layer,'String', text_slider); 

            [X,Z] = meshgrid(x*x_ab,z_start + z*z_ab);
            Figure = surf(X,Z,h_mat_xz(:,:,y_layer));
            
            text_title = strcat(title_label,' [ Y = ',' ',num2str((get(handles.layer_slider,'value')-(y_anzahl+1)/2)*y_ab),' mm ]');
            title(text_title)
            xlabel('x direction [mm]') 
            ylabel('z direction [mm]') 
            set(handles.axes1,'Xlim',x_limit,'Ylim',z_limit,'XTick',x_tick,'YTick',z_tick);
            view(2) % view 2D surf plot

            if interpolation == 1 
                shading interp
            else
                shading flat
            end
            
            colormap(jet)
            caxis(colorrange)
            
            %% Plot 3D figure 
            axes(handles.axes2);
            [X,Y,Z] = meshgrid(x*x_ab,y*y_ab,z_start + z*z_ab);
            slice(X,Y,Z,h_mat_xy,[],((get(handles.layer_slider,'value')-(y_anzahl+1)/2)*y_ab),[])
            set(handles.axes2,'Xlim',x_limit,'Ylim',y_limit,'Zlim',z_limit,'XTick',x_tick,'YTick',y_tick,'ZTick',z_tick)
            
            title(title_label);
            xlabel('x direction [mm]','Rotation',17) 
            ylabel('y direction [mm]','Rotation',-28) 
            zlabel('z direction [mm]') 
            
            colormap(jet)
            ylabel(colorbar, Colorbar_label)
            caxis(colorrange)
            
            if interpolation == 1 
                shading interp
            else
                shading flat
            end



        case 4

            set(handles.layer_slider,'Enable','on')
            textDebug = sprintf('');
            set(handles.debug, 'String', textDebug); 

            set(handles.layer_slider,'Enable','on')
            set(handles.layer_slider, 'min', x_limit_sb(1,1));
            set(handles.layer_slider, 'max', x_limit_sb(1,2));
                             
            set(handles.layer_slider, 'SliderStep',[1/(x_anzahl-1),1/(x_anzahl-1)]);           
            
            %% Plot 2d figure
            axes(handles.axes1);
            % set (gcf, 'WindowButtonMotionFcn', @mouseMove);
            textDebug = sprintf('');
            set(handles.debug, 'String', textDebug); 

            value = get(hObject,'Value');                  
            value = round(value);                                 
            set(hObject, 'Value', value);  
            x_layer =  get(handles.layer_slider,'value');                                            
            text_slider = sprintf('[ X = %d mm ]',(get(handles.layer_slider,'value')-(x_anzahl+1)/2)*x_ab);
            set(handles.text_layer,'String', text_slider); 
            [Y,Z] = meshgrid(y*y_ab,z_start + z*z_ab);
            Figure = surf(Y,Z,h_mat_yz(:,:,x_layer));
            text_title = strcat(title_label,' [ X = ',' ',num2str((get(handles.layer_slider,'value')-(x_anzahl+1)/2)*x_ab),' mm ]');
            title(text_title);
            xlabel('y direction [mm]') 
            ylabel('z direction [mm]')
            set(handles.axes1,'Xlim',y_limit,'Ylim',z_limit,'XTick',y_tick,'YTick',z_tick);
            view(2) % view 2D surf plot     

            if interpolation == 1 
                shading interp
            else
                shading flat
            end
            
            colormap(jet)
            caxis(colorrange)
            
            %% Plot 3D figure 
            axes(handles.axes2);
            [X,Y,Z] = meshgrid(x*x_ab,y*y_ab,z_start + z*z_ab);
            slice(X,Y,Z,h_mat_xy,((get(handles.layer_slider,'value')-(x_anzahl+1)/2)*x_ab),[],[])
            set(handles.axes2,'Xlim',x_limit,'Ylim',y_limit,'Zlim',z_limit,'XTick',x_tick,'YTick',y_tick,'ZTick',z_tick)
            
            title(title_label);
            xlabel('x direction [mm]','Rotation',17) 
            ylabel('y direction [mm]','Rotation',-28) 
            zlabel('z direction [mm]') 
            
            colormap(jet)
            ylabel(colorbar, Colorbar_label)
            caxis(colorrange)
            
            if interpolation == 1 
                shading interp
            else
                shading flat
            end
    end
else 
    msgbox('Please import the files and calculate the results first.');
    return
end










% --- Executes on selection change in data_type.
function data_type_Callback(hObject, eventdata, handles)
% hObject    handle to data_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns data_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from data_type

global flag
if flag == 2
    z_start =  str2double(get(handles.start_z ,'String'));
    x_anzahl = str2double(get(handles.x_anzahl ,'String'));
    y_anzahl  = str2double(get(handles.y_anzahl ,'String'));
    z_anzahl  = str2double(get(handles.z_anzahl ,'String'));

    x_ab = str2double(get(handles.x_abstand ,'String'));
    y_ab = str2double(get(handles.y_abstand ,'String'));
    z_ab = str2double(get(handles.z_abstand ,'String'));

    x_limit_sb = [1 x_anzahl];
    y_limit_sb = [1 y_anzahl];
    z_limit_sb = [1 z_anzahl];

    x_limit = [-(x_anzahl-1)*x_ab/2 (x_anzahl-1)*x_ab/2];
    y_limit = [-(y_anzahl-1)*y_ab/2 (y_anzahl-1)*y_ab/2];
    z_limit = [z_start z_start + (z_anzahl-1)*z_ab];
    
    x_tick = linspace(-(x_anzahl-1)*x_ab/2,(x_anzahl-1)*x_ab/2,4);
    y_tick = linspace(-(y_anzahl-1)*y_ab/2,(y_anzahl-1)*y_ab/2,4);
    z_tick = linspace(z_start, z_start + (z_anzahl-1)*z_ab,3);

    h_mat_xy = evalin('base','H_mat_xy');
    h_mat_xz = evalin('base','H_mat_xz');
    h_mat_yz = evalin('base','H_mat_yz');
    h_mat_xy_m = evalin('base','H_mat_xy_m');
    h_mat_xz_m = evalin('base','H_mat_xz_m');
    h_mat_yz_m = evalin('base','H_mat_yz_m');
    h_mat_xy_std = evalin('base','H_mat_xy_std');
    h_mat_xz_std = evalin('base','H_mat_xz_std');
    h_mat_yz_std = evalin('base','H_mat_yz_std');
    h_mat_xy_tb = evalin('base','H_mat_xy_tb');
    h_mat_xz_tb = evalin('base','H_mat_xz_tb');
    h_mat_yz_tb = evalin('base','H_mat_yz_tb');
    vel_rms = evalin('base','Vel_rms');

    interpolation =  get(handles.interpolation,'Value');
    plane = get(handles.plane,'Value');
    data_type = get(handles.data_type,'Value');
    
    x = -(x_anzahl-1)/2:(x_anzahl-1)/2;
    y = -(y_anzahl-1)/2:(y_anzahl-1)/2;
    z = 0 :(z_anzahl-1);  
    
     switch data_type
        case 1
             colorrange = [min(vel_rms(:,1)) max(vel_rms(:,1))];
             Colorbar_label = sprintf('Root Mean Square Velocity [m/s]');
             title_label =  sprintf('Root Mean Square Velocity ');
  
        case 2
             colorrange = [min(vel_rms(:,2)) max(vel_rms(:,2))];
             h_mat_xy = h_mat_xy_m;
             h_mat_xz = h_mat_xz_m;
             h_mat_yz = h_mat_yz_m;
             Colorbar_label = sprintf('Mean Velocity [m/s]');
             title_label =  sprintf('Mean Velocity ');
        case 3
             colorrange = [min(vel_rms(:,3)) max(vel_rms(:,3))];
             h_mat_xy = h_mat_xy_std;
             h_mat_xz = h_mat_xz_std;
             h_mat_yz = h_mat_yz_std;
             Colorbar_label = sprintf('Standard Deviation [m/s]');
             title_label =  sprintf('Standard Deviation ');
        case 4
             colorrange = [min(vel_rms(:,4)) max(vel_rms(:,4))];
             h_mat_xy = h_mat_xy_tb;
             h_mat_xz = h_mat_xz_tb;
             h_mat_yz = h_mat_yz_tb;
             Colorbar_label = sprintf('Turbulence Grad ');
             title_label =  sprintf('Turbulence Grad ');
    
    end

    switch plane

        case 1
            
            textDebug = sprintf('Please select a plane');
            set(handles.debug, 'String', textDebug); 
            set(handles.layer_slider,'Enable','off')
            set(handles.text_layer,'String', '');

        case 2
            textDebug = sprintf('');
            set(handles.debug, 'String', textDebug); 

            set(handles.layer_slider,'Enable','on')
            set(handles.layer_slider, 'min', z_limit_sb(1,1));
            set(handles.layer_slider, 'max', z_limit_sb(1,2));
                              
            set(handles.layer_slider, 'SliderStep',[1/(z_anzahl-1),1/(z_anzahl-1)]);           
            
            %% Plot 2D figure  
            axes(handles.axes1);
            textDebug = sprintf('');
            set(handles.debug, 'String', textDebug); 
            value = get(hObject,'Value');                  
            value = round(value);                                 
            set(hObject, 'Value', value);  
            z_layer =  get(handles.layer_slider,'value');                                            
            text_slider = sprintf('[ Z = %d mm ]',z_start + (z_layer-1)*z_ab);
            set(handles.text_layer,'String', text_slider); 
            [X,Y] = meshgrid(x*x_ab,y*y_ab);
            Figure = surf(X,Y,h_mat_xy(:,:,z_layer));
            text_title = strcat(title_label,' [ Z = ',' ',num2str(z_start + (z_layer-1)*z_ab),' mm ]');
            title(text_title);
            xlabel('x direction [mm]') 
            ylabel('y direction [mm]') 
            set(handles.axes1,'Xlim',x_limit,'Ylim',y_limit,'XTick',x_tick,'YTick',y_tick);
            view(2) % view 2D surf plot

            if interpolation == 1 
                shading interp
            else
                shading flat
            end
            
            colormap(jet)
            caxis(colorrange)
            

            %% Plot 3D figure 
            axes(handles.axes2);
            [X,Y,Z] = meshgrid(x*x_ab,y*y_ab,z_start + z*z_ab);
            slice(X,Y,Z,h_mat_xy,[],[],z_start + (z_layer-1)*z_ab)
            set(handles.axes2,'Xlim',x_limit,'Ylim',y_limit,'Zlim',z_limit,'XColor','r','YColor','[0.2 0 0]','ZColor','b','XTick',x_tick,'YTick',y_tick,'ZTick',z_tick)
            
            title(title_label);
            xlabel('x direction [mm]','Rotation',17) 
            ylabel('y direction [mm]','Rotation',-28) 
            zlabel('z direction [mm]') 
            
            colormap(jet)
            ylabel(colorbar, Colorbar_label)
            caxis(colorrange)
            
            if interpolation == 1 
                shading interp
            else
                shading flat
            end
            


        case 3

            set(handles.layer_slider,'Enable','on')
            textDebug = sprintf('');
            set(handles.debug, 'String', textDebug); 

            set(handles.layer_slider,'Enable','on')
            set(handles.layer_slider, 'min', y_limit_sb(1,1));
            set(handles.layer_slider, 'max', y_limit_sb(1,2));
                                
            set(handles.layer_slider, 'SliderStep',[1/(y_anzahl-1),1/(y_anzahl-1)]);           
            
            %% Plot 2D figure
            axes(handles.axes1);
            % set (gcf, 'WindowButtonMotionFcn', @mouseMove);
            textDebug = sprintf('');
            set(handles.debug, 'String', textDebug); 

            value = get(hObject,'Value');                  
            value = round(value);                                 
            set(hObject, 'Value', value);  
            y_layer =  get(handles.layer_slider,'value');                                            
            text_slider = sprintf('[ Y = %d mm ]',(get(handles.layer_slider,'value')-(y_anzahl+1)/2)*y_ab);
            set(handles.text_layer,'String', text_slider); 

            [X,Z] = meshgrid(x*x_ab,z_start + z*z_ab);
            Figure = surf(X,Z,h_mat_xz(:,:,y_layer));
            
            text_title = strcat(title_label,' [ Y = ',' ',num2str((get(handles.layer_slider,'value')-(y_anzahl+1)/2)*y_ab),' mm ]');
            title(text_title)
            xlabel('x direction [mm]') 
            ylabel('z direction [mm]') 
            set(handles.axes1,'Xlim',x_limit,'Ylim',z_limit,'XTick',x_tick,'YTick',z_tick);
            view(2) % view 2D surf plot

            if interpolation == 1 
                shading interp
            else
                shading flat
            end
            
            colormap(jet)
            caxis(colorrange)
            
            %% Plot 3D figure 
            axes(handles.axes2);
            [X,Y,Z] = meshgrid(x*x_ab,y*y_ab,z_start + z*z_ab);
            slice(X,Y,Z,h_mat_xy,[],((get(handles.layer_slider,'value')-(y_anzahl+1)/2)*y_ab),[])
            set(handles.axes2,'Xlim',x_limit,'Ylim',y_limit,'Zlim',z_limit,'XTick',x_tick,'YTick',y_tick,'ZTick',z_tick)
            
            title(title_label);
            xlabel('x direction [mm]','Rotation',17) 
            ylabel('y direction [mm]','Rotation',-28) 
            zlabel('z direction [mm]') 
            
            colormap(jet)
            ylabel(colorbar, Colorbar_label)
            caxis(colorrange)
            
            if interpolation == 1 
                shading interp
            else
                shading flat
            end



        case 4

            set(handles.layer_slider,'Enable','on')
            textDebug = sprintf('');
            set(handles.debug, 'String', textDebug); 

            set(handles.layer_slider,'Enable','on')
            set(handles.layer_slider, 'min', x_limit_sb(1,1));
            set(handles.layer_slider, 'max', x_limit_sb(1,2));
                               
            set(handles.layer_slider, 'SliderStep',[1/(x_anzahl-1),1/(x_anzahl-1)]);           
            
            %% Plot 2d figure
            axes(handles.axes1);
            % set (gcf, 'WindowButtonMotionFcn', @mouseMove);
            textDebug = sprintf('');
            set(handles.debug, 'String', textDebug); 

            value = get(hObject,'Value');                  
            value = round(value);                                 
            set(hObject, 'Value', value);  
            x_layer =  get(handles.layer_slider,'value');                                            
            text_slider = sprintf('[ X = %d mm ]',(get(handles.layer_slider,'value')-(x_anzahl+1)/2)*x_ab);
            set(handles.text_layer,'String', text_slider); 
            [Y,Z] = meshgrid(y*y_ab,z_start + z*z_ab);
            Figure = surf(Y,Z,h_mat_yz(:,:,x_layer));
            text_title = strcat(title_label,' [ X = ',' ',num2str((get(handles.layer_slider,'value')-(x_anzahl+1)/2)*x_ab),' mm ]');
            title(text_title);
            xlabel('y direction [mm]') 
            ylabel('z direction [mm]')
            set(handles.axes1,'Xlim',y_limit,'Ylim',z_limit,'XTick',y_tick,'YTick',z_tick);
            view(2) % view 2D surf plot     

            if interpolation == 1 
                shading interp
            else
                shading flat
            end
            
            colormap(jet)
            caxis(colorrange)
            
            %% Plot 3D figure 
            axes(handles.axes2);
            [X,Y,Z] = meshgrid(x*x_ab,y*y_ab,z_start + z*z_ab);
            slice(X,Y,Z,h_mat_xy,((get(handles.layer_slider,'value')-(x_anzahl+1)/2)*x_ab),[],[])
            set(handles.axes2,'Xlim',x_limit,'Ylim',y_limit,'Zlim',z_limit,'XTick',x_tick,'YTick',y_tick,'ZTick',z_tick)
            
            title(title_label);
            xlabel('x direction [mm]','Rotation',17) 
            ylabel('y direction [mm]','Rotation',-28) 
            zlabel('z direction [mm]') 
            
            colormap(jet)
            ylabel(colorbar, Colorbar_label)
            caxis(colorrange)
            
            if interpolation == 1 
                shading interp
            else
                shading flat
            end
    end
else 
    msgbox('Please import the files and calculate the results first.');
    return
end











% --- Executes during object creation, after setting all properties.
function data_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function interpolation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to interpolation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function layer_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to layer_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
   
% --- Executes on selection change in Z Direction.
function eben_Callback(hObject, eventdata, handles)
% hObject    handle to eben (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns eben contents as cell array
%        contents{get(hObject,'Value')} returns selected item from eben


% --- Executes during object creation, after setting all properties.
function eben_CreateFcn(hObject, eventdata, handles)
% hObject    handle to eben (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in exportieren.
function exportieren_Callback(hObject, eventdata, handles)
% hObject    handle to exportieren (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes1 = handles.axes1;
saveFig = figure('visible','off');
saveData = copyobj(axes1,saveFig);
set(saveData,'Units','default','Position','default');
colormap(jet);
ylabel(colorbar,'Velocity [m/s]');
set(gcf,'Name','Preview of the figure','NumberTitle','off');

[filename, pathname, ~] = uiputfile( {'*.fig';'*.png';'*.emf';'*.jpg';'*.tif';},'save the figure','untitled');
     if filename==0
        % user pressed cancel
        close(saveFig);
        delete(saveData);
        return

     end
     
         set(saveFig,'Visible','On')
         savedata = fullfile(pathname, filename);
         saveas(saveData,savedata);
%          close(saveFig);
%          delete(saveData);
    
    
 
 
 
 
 % --- Executes on button press in export_3d.
function export_3d_Callback(hObject, eventdata, handles)
% hObject    handle to export_3d (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes1 = handles.axes2;
saveFig = figure('visible','off');
saveData = copyobj(axes1,saveFig);
set(saveData,'Units','default','Position','default');
colormap(jet);
ylabel(colorbar,'Velocity [m/s]');
set(gcf,'Name','Preview of the figure','NumberTitle','off');

[filename, pathname, ~] = uiputfile( {'*.fig';'*.png';'*.emf';'*.jpg';'*.tif';},'save the figure','untitled_3D');
     if filename==0
        % user pressed cancel
        close(saveFig);
        delete(saveData);
        return

     end
         set(saveFig,'Visible','On')
         savedata = fullfile(pathname, filename);
         saveas(saveData,savedata);
%          close(saveFig);
%          delete(saveData);

    
     
     
     % --- Executes on button press in export_data.
function export_data_Callback(hObject, eventdata, handles)
% hObject    handle to export_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% cell_coor =  mat2cell(coordinate,ones(1,size(coordinate,1)),ones(1,5));
% cell_result = [{'x [mm]','y [mm]','z [mm]','vStr,mittel [m/s]','Standardabweichung'};cell_coor];
global flag
if flag == 2
vel_rms = evalin('base','Vel_rms');
data =  result(handles, vel_rms);
uisave({'data'},'')
else 
    msgbox('Please import the files and calculate the results first.');
    return
end




function Temperatur_Callback(hObject, eventdata, handles)
% hObject    handle to Temperatur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Temperatur as text
%        str2double(get(hObject,'String')) returns contents of Temperatur as a double


% --- Executes during object creation, after setting all properties.
function Temperatur_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Temperatur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






function x_abstand_Callback(hObject, eventdata, handles)
% hObject    handle to x_abstand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x_abstand as text
%        str2double(get(hObject,'String')) returns contents of x_abstand as a double


% --- Executes during object creation, after setting all properties.
function x_abstand_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_abstand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y_abstand_Callback(hObject, eventdata, handles)
% hObject    handle to y_abstand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y_abstand as text
%        str2double(get(hObject,'String')) returns contents of y_abstand as a double


% --- Executes during object creation, after setting all properties.
function y_abstand_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_abstand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function z_abstand_Callback(hObject, eventdata, handles)
% hObject    handle to z_abstand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of z_abstand as text
%        str2double(get(hObject,'String')) returns contents of z_abstand as a double


% --- Executes during object creation, after setting all properties.
function z_abstand_CreateFcn(hObject, eventdata, handles)
% hObject    handle to z_abstand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function x_anzahl_Callback(hObject, eventdata, handles)
% hObject    handle to x_anzahl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of x_anzahl as text
%        str2double(get(hObject,'String')) returns contents of x_anzahl as a double


% --- Executes during object creation, after setting all properties.
function x_anzahl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to x_anzahl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function z_anzahl_Callback(hObject, eventdata, handles)
% hObject    handle to z_anzahl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of z_anzahl as text
%        str2double(get(hObject,'String')) returns contents of z_anzahl as a double


% --- Executes during object creation, after setting all properties.
function z_anzahl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to z_anzahl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function y_anzahl_Callback(hObject, eventdata, handles)
% hObject    handle to y_anzahl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of y_anzahl as text
%        str2double(get(hObject,'String')) returns contents of y_anzahl as a double


% --- Executes during object creation, after setting all properties.
function y_anzahl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to y_anzahl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function debug_CreateFcn(hObject, eventdata, handles)
% hObject    handle to debug (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function plane_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function layer_Callback(hObject, eventdata, handles)
% hObject    handle to layer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of layer as text
%        str2double(get(hObject,'String')) returns contents of layer as a double


% --- Executes during object creation, after setting all properties.
function layer_CreateFcn(hObject, eventdata, handles)
% hObject    handle to layer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in start_pot.
function start_pot_Callback(hObject, eventdata, handles)
% hObject    handle to start_pot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns start_pot contents as cell array
%        contents{get(hObject,'Value')} returns selected item from start_pot


% --- Executes during object creation, after setting all properties.
function start_pot_CreateFcn(hObject, eventdata, handles)
% hObject    handle to start_pot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in start_direct.


function start_direct_Callback(hObject, eventdata, handles)
% hObject    handle to start_direct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns start_direct contents as cell array
%        contents{get(hObject,'Value')} returns selected item from start_direct


% --- Executes during object creation, after setting all properties.
function start_direct_CreateFcn(hObject, eventdata, handles)
% hObject    handle to start_direct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function text_layer_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_layer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function start_x_Callback(hObject, eventdata, handles)
% hObject    handle to start_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of start_x as text
%        str2double(get(hObject,'String')) returns contents of start_x as a double


% --- Executes during object creation, after setting all properties.
function start_x_CreateFcn(hObject, eventdata, handles)
% hObject    handle to start_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function start_y_Callback(hObject, eventdata, handles)
% hObject    handle to start_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of start_y as text
%        str2double(get(hObject,'String')) returns contents of start_y as a double


% --- Executes during object creation, after setting all properties.
function start_y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to start_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function start_z_Callback(hObject, eventdata, handles)
% hObject    handle to start_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of start_z as text
%        str2double(get(hObject,'String')) returns contents of start_z as a double


% --- Executes during object creation, after setting all properties.
function start_z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to start_z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a5_Callback(hObject, eventdata, handles)
% hObject    handle to a5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a5 as text
%        str2double(get(hObject,'String')) returns contents of a5 as a double


% --- Executes during object creation, after setting all properties.
function a5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a2_Callback(hObject, eventdata, handles)
% hObject    handle to a2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a2 as text
%        str2double(get(hObject,'String')) returns contents of a2 as a double


% --- Executes during object creation, after setting all properties.
function a2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a4_Callback(hObject, eventdata, handles)
% hObject    handle to a4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a4 as text
%        str2double(get(hObject,'String')) returns contents of a4 as a double


% --- Executes during object creation, after setting all properties.
function a4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a1_Callback(hObject, eventdata, handles)
% hObject    handle to a1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a1 as text
%        str2double(get(hObject,'String')) returns contents of a1 as a double


% --- Executes during object creation, after setting all properties.
function a1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a3_Callback(hObject, eventdata, handles)
% hObject    handle to a3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a3 as text
%        str2double(get(hObject,'String')) returns contents of a3 as a double


% --- Executes during object creation, after setting all properties.
function a3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function a0_Callback(hObject, eventdata, handles)
% hObject    handle to a0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of a0 as text
%        str2double(get(hObject,'String')) returns contents of a0 as a double


% --- Executes during object creation, after setting all properties.
function a0_CreateFcn(hObject, eventdata, handles)
% hObject    handle to a0 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function sample_rate_Callback(hObject, eventdata, handles)
% hObject    handle to sample_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sample_rate as text
%        str2double(get(hObject,'String')) returns contents of sample_rate as a double


% --- Executes during object creation, after setting all properties.
function sample_rate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sample_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
