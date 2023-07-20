function exportVTK(pathToFile,keyExport,element,mesh,U,V,STRESS,T,HMH,eps_eq,problem)

etyp = element.type;
X = mesh.X;
IEN = mesh.IEN;

nnp = size(X,2);
nsd = size(X,1);
if(nsd == 2)
    X = [X; zeros(1,nnp)];
end

nel = size(IEN,2);
nen = size(IEN,1);

if(~isempty(U))
    nnd = size(U,1);
    if(nnd == 2)
        U = [U; zeros(1,nnp)];
    end
end

if(~isempty(V))
    nnd = size(V,1);
    if(nnd == 2)
        V = [V; zeros(1,nnp)];
    end
end

if(~isempty(STRESS))
    if(size(STRESS,1) == 3)        
        STRESS = [STRESS(1,:); STRESS(2,:); STRESS(3,:); zeros(1,nnp); zeros(1,nnp); zeros(1,nnp)];
    end
    if(size(STRESS,1) == 4)        
        STRESS = [STRESS(1,:); STRESS(2,:); STRESS(3,:); STRESS(4,:); zeros(1,nnp); zeros(1,nnp)];
    end
end

%%
fout = fopen(pathToFile,'w');
fprintf(fout,'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n');
fprintf(fout,'  <UnstructuredGrid>\n');
fprintf(fout,'    <Piece NumberOfPoints="%i" NumberOfCells="%i">\n',nnp,nel);

fprintf(fout,'	  <Points>\n');
fprintf(fout,'        <DataArray type="Float64" NumberOfComponents="3" format="ascii">\n');
fprintf(fout,'          %f %f %f\n',X);
fprintf(fout,'        </DataArray>\n');
fprintf(fout,'	  </Points>\n');

fprintf(fout,'      <Cells>\n');
fprintf(fout,'		  <DataArray type="Int32" Name="connectivity" format="ascii">\n');
switch etyp
    case 'p1quadD'
        VTK_CODE = 9;
        fprintf(fout,'          %d %d %d %d\n',(IEN(1:nen,:)-1));
     case 'p1shellD'
        VTK_CODE = 9;
        fprintf(fout,'          %d %d %d %d\n',(IEN(1:nen,:)-1));       
    case 'p1quadT'
        VTK_CODE = 9;
        fprintf(fout,'          %d %d %d %d\n',(IEN(1:nen,:)-1));            
    case 'p1hexT'
        VTK_CODE = 12;
        fprintf(fout,'          %d %d %d %d %d %d %d %d\n',(IEN(1:nen,:)-1));        
    case 'p1hexD'
        VTK_CODE = 12;
        fprintf(fout,'          %d %d %d %d %d %d %d %d\n',(IEN(1:nen,:)-1));        
    case 'p2quad'
        VTK_CODE = 23;
        fprintf(fout,'          %d %d %d %d %d %d %d %d\n',(IEN(1:nen,:)-1));
    case 'p2quadD'
        VTK_CODE = 23;
        fprintf(fout,'          %d %d %d %d %d %d %d %d\n',(IEN(1:nen,:)-1));
    case 'p2hexD'
        VTK_CODE = 25;
        fprintf(fout,'          %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n',(IEN([1:12 17:20 13:16],:)-1));
    case 'p2hexT'
        VTK_CODE = 25;
        fprintf(fout,'          %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n',(IEN([1:12 17:20 13:16],:)-1));        
    otherwise
        warning('Unexpected element type.');
end

fprintf(fout,'        </DataArray>\n');
fprintf(fout,'        <DataArray type="Int32" Name="offsets" format="ascii">\n');
fprintf(fout,'          %d\n',(nen:nen:nen*nel));
fprintf(fout,'        </DataArray>\n');
fprintf(fout,'        <DataArray type="UInt8" Name="types" format="ascii">\n');
fprintf(fout,'          %d\n',ones(1,nel)*VTK_CODE);
fprintf(fout,'        </DataArray>\n');
fprintf(fout,'      </Cells>\n');

fprintf(fout,'      <PointData Scalars="scalars">\n');
if(keyExport.displacement)
    fprintf(fout,'        <DataArray type="Float64" Name="displacement" NumberOfComponents="3" Format="ascii">\n');
    fprintf(fout,'          %.8e %.8e %.8e\n',full(U));
    fprintf(fout,'        </DataArray>\n');
end

if(keyExport.velocity)
    fprintf(fout,'        <DataArray type="Float64" Name="velocity" NumberOfComponents="3" Format="ascii">\n');
    fprintf(fout,'          %.8e %.8e %.8e\n',full(V));
    fprintf(fout,'        </DataArray>\n');
end

if(keyExport.temperature)
    fprintf(fout,'        <DataArray type="Float64" Name="temperature" NumberOfComponents="1" Format="ascii">\n');
    fprintf(fout,'          %.8e\n',full(T)');
    fprintf(fout,'        </DataArray>\n');
end

%fprintf(fout,'        <DataArray type="Float64" Name="HMH" NumberOfComponents="1" Format="ascii">\n');
%fprintf(fout,'          %.8f\n',HMH');
%fprintf(fout,'        </DataArray>\n');

if(keyExport.stress)
	fprintf(fout,'        <DataArray type="Float64" Name="stress" NumberOfComponents="6" Format="ascii">\n');
	fprintf(fout,'          %.8f %.8f %.8f %.8f %.8f %.8f\n',STRESS);
	fprintf(fout,'        </DataArray>\n');
end

if(keyExport.HMH)
    fprintf(fout,'        <DataArray type="Float64" Name="HMH" NumberOfComponents="1" Format="ascii">\n');
    fprintf(fout,'          %.8e\n',full(HMH)');
    fprintf(fout,'        </DataArray>\n');
end

if(keyExport.eps_eq)
    fprintf(fout,'        <DataArray type="Float64" Name="eps_eqv" NumberOfComponents="1" Format="ascii">\n');
    fprintf(fout,'          %.8e\n',full(eps_eq)');
    fprintf(fout,'        </DataArray>\n');
end

fprintf(fout,'      </PointData>\n');

if(keyExport.partitioning)

partition_ID = load(strcat(problem.dir,problem.meshName,'.epart.',problem.nuberOfPartitions));

fprintf(fout,'      <CellData Scalars="scalars">\n');
fprintf(fout,'        <DataArray type="Int32" Name="partition" NumberOfComponents="1" Format="ascii">\n');
fprintf(fout,'          %d\n',full(partition_ID)');
fprintf(fout,'        </DataArray>\n');
fprintf(fout,'      </CellData>\n');
end

fprintf(fout,'    </Piece>\n');
fprintf(fout,'  </UnstructuredGrid>\n');
fprintf(fout,'</VTKFile>\n');

fclose(fout);
