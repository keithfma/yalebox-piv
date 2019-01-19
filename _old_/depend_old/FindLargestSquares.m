function S = FindLargestSquares(I)
%FindLargestSquares - finds largest sqare regions with all points set to 1.
%input:  I - B/W boolean matrix
%output: S - for each pixel I(r,c) return size of the largest all-white square with its upper -left corner at I(r,c)  

% From Inscribed_Rectangle by Jaroslaw Tuszynski 08 Jul 2010 (Updated 06 Jan 2012) 
% http://www.mathworks.com/matlabcentral/fileexchange/28155-inscribed-rectangle
%
% LICENSE:
% Copyright (c) 2010, Jaroslaw Tuszynski
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

[nr, nc] = size(I);
S = double(I>0);
for r=(nr-1):-1:1
  for c=(nc-1):-1:1
    if (S(r,c))
      a = S(r  ,c+1);
      b = S(r+1,c  );
      d = S(r+1,c+1);
      S(r,c) = min([a b d]) + 1;
    end
  end  
end  