% This is for .xyz coordinates file format with only sphere for representing atoms
clc
clear all
%% Read lines of the text file
% Name of the .xyz file
A = readlines("ACF.dat");
%% Inputs
% Reading the number of atoms from ACF.dat file. Comfirm that the 5th line
% from end of page is the description of last atom
NAtoms_str = A(end-5);
NAtoms_dou = ReadLine(NAtoms_str, 1);
NAtoms = NAtoms_dou(1); % Number of atoms in the ACF.dat file
% Initializing the variables taken from POSCAR file
AtomIndex = zeros(NAtoms,1); % Index of the atoms
Coord = zeros(NAtoms,3); % X,Y,Z coordinates of the atoms
Charge = zeros(NAtoms,1); % The total electrons obtained from ACF.dat file
AtomColor = zeros(NAtoms, 3);
ValElec = zeros(NAtoms, 1);
AtomSize = ones(NAtoms, 1);
%% Reading the POSCAR file if it exists
if (isfile("POSCAR"))
    B = readlines("POSCAR");
    % Getting the involved atoms in Line 1
    Elem = ReadLine(B(1), 0);
    NumAtomElem = ReadLine(B(6), 1);
    CurrInd = 1;
    % Specifying the color, number of valence electrons and current index
    % for the atoms. If your atom is not included in the list, those values
    % will be initialized to 0 and you would not be able to differentiate the atoms.
    for i=1:length(Elem)
        if Elem(i) == "C"
            AtomColor(CurrInd:CurrInd+NumAtomElem(i)-1, :) = repmat([0.5 0.5 0.5], NumAtomElem(i), 1);
            ValElec(CurrInd:CurrInd+NumAtomElem(i)-1, :) = ones(NumAtomElem(i),1)*4;
            AtomSize(CurrInd:CurrInd+NumAtomElem(i)-1, :) = ones(NumAtomElem(i),1)*0.5;
            CurrInd = CurrInd+NumAtomElem(i);
        elseif Elem(i) == "F"
            AtomColor(CurrInd:CurrInd+NumAtomElem(i)-1, :) = repmat([179/255 1 1], NumAtomElem(i), 1);
            ValElec(CurrInd:CurrInd+NumAtomElem(i)-1, :) = ones(NumAtomElem(i),1)*7;
            AtomSize(CurrInd:CurrInd+NumAtomElem(i)-1, :) = ones(NumAtomElem(i),1)*0.5;
            CurrInd = CurrInd+NumAtomElem(i);
        elseif Elem(i) == "N"
            AtomColor(CurrInd:CurrInd+NumAtomElem(i)-1, :) = repmat([0 0 1], NumAtomElem(i), 1);
            ValElec(CurrInd:CurrInd+NumAtomElem(i)-1, :) = ones(NumAtomElem(i),1)*5;
            AtomSize(CurrInd:CurrInd+NumAtomElem(i)-1, :) = ones(NumAtomElem(i),1)*0.5;
            CurrInd = CurrInd+NumAtomElem(i);
        elseif Elem(i) == "O"
            AtomColor(CurrInd:CurrInd+NumAtomElem(i)-1, :) = repmat([1 0. 0.], NumAtomElem(i), 1);
            ValElec(CurrInd:CurrInd+NumAtomElem(i)-1, :) = ones(NumAtomElem(i),1)*6;
            AtomSize(CurrInd:CurrInd+NumAtomElem(i)-1, :) = ones(NumAtomElem(i),1)*0.5;
            CurrInd = CurrInd+NumAtomElem(i);
        elseif Elem(i) == "H"
            AtomColor(CurrInd:CurrInd+NumAtomElem(i)-1, :) = repmat([0.9 0.9 0.9], NumAtomElem(i), 1);
            ValElec(CurrInd:CurrInd+NumAtomElem(i)-1, :) = ones(NumAtomElem(i),1)*1;
            AtomSize(CurrInd:CurrInd+NumAtomElem(i)-1, :) = ones(NumAtomElem(i),1)*0.2;
            CurrInd = CurrInd+NumAtomElem(i);
        elseif Elem(i) == "Pt"
            AtomColor(CurrInd:CurrInd+NumAtomElem(i)-1, :) = repmat([24/255 91/255 145/255], NumAtomElem(i), 1);
            ValElec(CurrInd:CurrInd+NumAtomElem(i)-1, :) = ones(NumAtomElem(i),1)*10;
            AtomSize(CurrInd:CurrInd+NumAtomElem(i)-1, :) = ones(NumAtomElem(i),1)*1;
            CurrInd = CurrInd+NumAtomElem(i);
        elseif Elem(i) == "Ru"
            AtomColor(CurrInd:CurrInd+NumAtomElem(i)-1, :) = repmat([35/255 142/255 151/255], NumAtomElem(i), 1);
            ValElec(CurrInd:CurrInd+NumAtomElem(i)-1, :) = ones(NumAtomElem(i),1)*8;
            AtomSize(CurrInd:CurrInd+NumAtomElem(i)-1, :) = ones(NumAtomElem(i),1)*1;
            CurrInd = CurrInd+NumAtomElem(i);
        elseif Elem(i) == "Fe"
            AtomColor(CurrInd:CurrInd+NumAtomElem(i)-1, :) = repmat([129/255 123/255 196/255], NumAtomElem(i), 1);
            ValElec(CurrInd:CurrInd+NumAtomElem(i)-1, :) = ones(NumAtomElem(i),1)*8;
            AtomSize(CurrInd:CurrInd+NumAtomElem(i)-1, :) = ones(NumAtomElem(i),1)*1;
            CurrInd = CurrInd+NumAtomElem(i);
        elseif Elem(i) == "Cu"
            AtomColor(CurrInd:CurrInd+NumAtomElem(i)-1, :) = repmat([1 123/255 98/255], NumAtomElem(i), 1);
            ValElec(CurrInd:CurrInd+NumAtomElem(i)-1, :) = ones(NumAtomElem(i),1)*11;
            AtomSize(CurrInd:CurrInd+NumAtomElem(i)-1, :) = ones(NumAtomElem(i),1)*1;
            CurrInd = CurrInd+NumAtomElem(i);
        elseif Elem(i) == "Ag"
            AtomColor(CurrInd:CurrInd+NumAtomElem(i)-1, :) = repmat([1 235/255 0], NumAtomElem(i), 1);
            ValElec(CurrInd:CurrInd+NumAtomElem(i)-1, :) = ones(NumAtomElem(i),1)*6;
            AtomSize(CurrInd:CurrInd+NumAtomElem(i)-1, :) = ones(NumAtomElem(i),1)*1;
            CurrInd = CurrInd+NumAtomElem(i);
        end
    end
end
%% Extract all the data from ACF.dat file and store in arrays
% The first two lines are not needed. Start at line 3.
for i=1:NAtoms
    dum2 = A(i+2); % The string line containing the atom index, coordinates and charges
    C = ReadLine(dum2, 1);
    AtomIndex(i) = C(1);
    Coord(i, 1:3) = C(2:4);
    Charge(i) = C(5);
end
%% Plotting metal and adsorbate atoms in a scatter plot
f = figure(1);
axis equal;
hScatter = scatter3(Coord(:,1), Coord(:,2), Coord(:,3), 600*AtomSize, AtomColor, 'filled');
% If you want to enable charges visibility only on hovering
dcm = datacursormode(gcf);
set(dcm, 'UpdateFcn', @(obj, event) myTooltip(event, string(ValElec-Charge)));
% If you want to enable all the charge labels visible at all times
% dx = 0.1;
% text(Coord(:,1)+dx, Coord(:,2)+dx, Coord(:,3)+dx, string(Charge));
%% Functions
% Get all the words in an array
function Coords = ReadLine(A, toDouble)
    Coords = [];
    NumbChar = strlength(A);
    i = 1;
    while i<=NumbChar
        dum1 = extract(A,i);
        word = "";
        EnteredLoop = 0;
        while dum1 ~= " "
            EnteredLoop = 1;
            if (i<NumbChar)
                word = word + dum1;
                i = i + 1;
                dum1 = extract(A,i);
            else
                word = word + dum1;
                break;
            end
        end
        if (EnteredLoop == 1)
            if (toDouble)
                Coords = [Coords str2double(word)];
            else
                Coords = [Coords word];
            end
        end
        i = i + 1;
    end
end
function output_txt = myTooltip(event, labels)
    % Get index of the selected data point
    idx = event.DataIndex;
    % Show coordinates and custom label
    output_txt = {[labels{idx}]};
end