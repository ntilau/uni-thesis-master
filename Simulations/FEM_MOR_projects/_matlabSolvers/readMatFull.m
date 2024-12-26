function Mat = readMatFull( filename )

fid = fopen( filename, 'r' );
if fid == -1
  error(strcat('Could not open file: ', filename));
end

% read '{' out of file
blockBegin = fscanf(fid, '%s', 1);
if ~strcmp(blockBegin, '{')
  error('File doesnt begin with {');
end

matType = fscanf( fid, '%s', 1 );

if matType(4)=='S'  % matrix is symmetric
  dim = fscanf( fid, '%i', 1 );
  Mat = zeros(dim,dim);
  if matType(2) == 'C'
    for col=1:dim
      for row=1:col
        cVal = readComplex(fid);
        Mat(row,col) = cVal;
        Mat(col,row) = cVal;
      end
    end
  elseif matType(2) == 'R'
    for col=1:dim
      for row=1:col
        val = fscanf(fid, '%f', 1);
        Mat(row,col) = val;
        Mat(col,row) = val;
      end
    end
  else
    error('Something is wrong with the matrix type');
  end
elseif matType(4)=='R'  % matrix is rectangular
  numRows = fscanf( fid, '%i', 1 );
  numCols = fscanf( fid, '%i', 1 );
  Mat = zeros(numRows,numCols);
  if matType(2) == 'C'
    for row=1:numRows
      for col=1:numCols
        val = readComplex(fid);
        Mat(row,col) = val;
      end
    end
    % error('Not yet implemented');
  elseif matType(2) == 'R'
    for row=1:numRows
      for col=1:numCols
        val = fscanf(fid, '%g', 1);
        Mat(row,col) = val;
      end
    end
  else
    error('Something is wrong with the matrix type'); 
  end
elseif matType(4)=='Q'  % matrix is quadratic
  dim = fscanf( fid, '%i', 1 );
  Mat = zeros(dim);
  if matType(2) == 'C'
    for row=1:dim
      for col=1:dim
        val = readComplex(fid);
        Mat(row,col) = val;
      end
    end
    % error('Not yet implemented');
  elseif matType(2) == 'R'
    for row=1:dim
      for col=1:dim
        val = fscanf(fid, '%g', 1);
        Mat(row,col) = val;
      end
    end
  else
    error('Something is wrong with the matrix type'); 
  end
else
  error('Not yet implemented');
end
fclose(fid);