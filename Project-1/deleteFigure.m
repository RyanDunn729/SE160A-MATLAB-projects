function result = deleteFigure(fname)
% DELETEFIGURE Delete all figures in an Excel file.
%   
%   Usage:
%   DELETEFIGURE()          will ask for an input file name to delete all figures in.
%   DELETEFIGURE(fname)     deletes all figures in fname file.
%   A = DELETEFIGURE()      will ask for an input file name to delete all figures in.
%   A = DELETEFIGURE(fname) deletes all figures in fname file.
%   
%   The optional output A will return 1 if successful, 0 if not.
%   
%   Example 1:
%   A = DELETEFIGURE();
%   Please input Excel file name.
%   ?> ExcelFileName.xlsx
%   
%   Example 2:
%   DELETEFIGURE('ExcelFileName.xlsx');

% -------------------------------------------------------------------------
%   UCSD - SE160A / 260A - Professor John B. Kosmatka
%   Written by: James Chen
%   Last revised: 2016-03-03
% -------------------------------------------------------------------------


if(nargin == 0); % if no input argument
    fprintf('Please input Excel file name.\n');
    fname = input('?>\ ','s');
end

file = which(fname); % get full path of file
if(isempty(file)); % if undefined file
    result = 0;
    fprintf('Undefined file.\n');
    return;
end

try
    Excel = actxserver('Excel.Application'); % start activex server
catch
    fprintf('Cannot start ActiveX server.\n');
    result = 0;
    return;
end

try
    Workbook = Excel.Workbooks.Open(file); % open file
    Worksheet = Workbook.sheets; % get list of sheets
    numSheets = Worksheet.Count; % number of sheets
    for i = 1:numSheets; % loop over all sheets
        Shape = Worksheet.Item(i).Shapes; % get shape data
        for j = Shape.Count:-1:1; % loop over all shapes within sheet
            Shape.Item(j).Delete; % delete current shape
        end
    end
    
    Workbook.Save; % save file
    Workbook.Close;
    Excel.Quit;
    result = 1; % successful result
catch
    fprintf('Error.\n'); % if something errored
    Excel.Quit;
    result = 0; % unsuccessful result
end


end