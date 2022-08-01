function result = createFigure(fname,sheetNum,figNum,TL,BR)
% CREATEFIGURE Output a plot figure to an Excel file with defined
% boundaries. 
%   
%   Usage:
%   CREATEFIGURE(fname,sheetNum,figNum,topLeft,bottomRight);
%   A = CREATEFIGURE(fname,sheetNum,figNum,topLeft,bottomRight);
%   
%   The optional output A will return 1 if successful, 0 if not.
%   fname        name of Excel file (e.g. 'MyExcelOutput.xlsx')
%                File can include relative or absolute path to file.
%   sheetNum     sheet number (e.g. 1, 2, or 3)
%                Sheet must already exist within the file.
%   figNum       figure handle to output
%                Figure must have been pre-iniliatized prior to calling
%                this function.
%   topLeft      top-left cell bound (e.g. 'B2')
%   bottomRight  bottom-right cell bound (e.g. 'H20')
%   
%   Example 1:
%   x = 0:0.1:10; y = sin(x);
%   figure(1); plot(x,y);
%   A = CREATEFIGURE('MyExcelOutput.xlsx', 1, 1, 'B2', 'H20');
%   
%   Example 2:
%   x = -10:0.1:10; y = x.^2;
%   figure(3); plot(x,y);
%   CREATEFIGURE('..\MyExcelOutput.xlsx', 1, 3, 'B2', 'H20');

% -------------------------------------------------------------------------
%   UCSD - SE160A / 260A - Professor John B. Kosmatka
%   Written by: James Chen
%   Last revised: 2016-03-03
% -------------------------------------------------------------------------


img = 'tempHolder.png'; % temporary image file placeholder

file = which(fname); % get full path of file

% Check if creation of a new file is necessary.
createNew = 0;
if(isempty(file)); % if undefined file
    fprintf('File not found. Creating new file.\n');
    createNew = 1; % then create new file
end

% Start ActiveX server (required for this function).
try
    Excel = actxserver('Excel.Application'); % start activex server
catch
    fprintf('Cannot start ActiveX server.\n');
    result = 0;
    return; % abort if failed
end

try
    if(~createNew); % if file already exists
        Workbook = Excel.Workbooks.Open(file); % open file
    else
        Workbook = invoke(Excel.Workbooks,'Add'); % create new file
        sheetNum = 1; % default to first sheet
    end
    
    Sheets = Workbook.Sheets; % get list of sheets
    Sheet = get(Sheets,'Item',sheetNum); % go to sheet specified by input
    Range = Sheet.Range(TL); % top left range
    Range2 = Sheet.Range(BR); % bottom right range
    
    % Geometrical calculations for bounding.
    StartRow = Range.Row;
    StartCol = Range.Column;
    EndRow = Range2.Row;
    EndCol = Range2.Column;
    
    if(EndRow <= StartRow || EndCol <= StartCol); % if negative or zero dimensions
        fprintf('Invalid range.\n');
        result = 0;
        return;
    end
    
    TotalHeight = 0; % height of image
    for i = StartRow:1:EndRow;
        tempAddress = get(get(Sheet,'Cells',i,StartCol),'Address');
        tempRange = Sheet.Range(tempAddress);
        TotalHeight = TotalHeight + tempRange.RowHeight;
    end
    
    TotalWidth = 0; % width of image
    for i = StartCol:1:EndCol;
        tempAddress = get(get(Sheet,'Cells',StartRow,i),'Address');
        tempRange = Sheet.Range(tempAddress);
        TotalWidth = TotalWidth + tempRange.ColumnWidth;
    end
    
    OffsetHeight = 0; % top bound
    for i = 1:1:StartRow-1;
        tempAddress = get(get(Sheet,'Cells',i,StartCol),'Address');
        tempRange = Sheet.Range(tempAddress);
        OffsetHeight = OffsetHeight + tempRange.RowHeight;
    end
    
    OffsetWidth = 0; % left bound
    for i = 1:1:StartCol-1;
        tempAddress = get(get(Sheet,'Cells',StartRow,i),'Address');
        tempRange = Sheet.Range(tempAddress);
        OffsetWidth = OffsetWidth + tempRange.ColumnWidth;
    end
    
    b = 5/3*0.72/0.21; % cell size ratio
    
    Left = OffsetWidth*b;
    Top = OffsetHeight;
    Width = TotalWidth*b;
    Height = TotalHeight;
    
    Shape = Sheet.Shapes;
    figure(figNum); % activate figure specified
    print('-dpng',img); % save to temporary image file
    Shape.AddPicture([pwd,filesep,img],0,1,Left,Top,Width,Height); % add image to file
    delete([pwd,filesep,img]); % delete temporary image file
    
    if(~createNew); % if file already exists
        Workbook.Save; % save existing file
    else % create new file with specified fname
        [path,name,ext] = fileparts(fname);
        if(isempty(path)); % if path was not specified
            invoke(Workbook,'SaveAs',[pwd,filesep,fname]); % save to current directory
        else % save to specified path
            curDir = pwd;
            cd(path);
            path = pwd;
            cd(curDir);
            invoke(Workbook,'SaveAs',[path,filesep,name,ext]);
        end
    end
    
    Workbook.Close;
    Excel.Quit;
    result = 1; % successful result
catch
    fprintf('Error.\n'); % if something errored
    Excel.Quit;
    result = 0; % unsuccessful result
end

end
