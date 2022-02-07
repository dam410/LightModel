% Copyright �2013. The Regents of the University of California (Regents).
% All Rights Reserved. Permission to use, copy, modify, and distribute
% this software and its documentation for educational, research, and
% not-for-profit purposes, without fee and without a signed licensing
% agreement, is hereby granted, provided that the above copyright notice,
% this paragraph and the following two paragraphs appear in all copies,
% modifications, and distributions. Contact The Office of Technology
% Licensing, UC Berkeley, 2150 Shattuck Avenue, Suite 510, Berkeley, CA
% 94720-1620, (510) 643-7201, for commercial licensing opportunities.
%
% Created by Jonathan T Barron and Jitendra Malik, Electrical Engineering
% and Computer Science, University of California, Berkeley.
%
% IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
% SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
% ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
% REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
% REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
% PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY,
% PROVIDED HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO
% PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.




f1 = [1, 2, 1;
      0, 0, 0;
    -1, -2, -1]/8;

f2 = [1, 0, -1;
      2, 0, -2;
      1, 0, -1]/8;


f12 = [1, 0, -1;
       0, 0,  0;
      -1, 0,  1]/4;


f11 =  [1,  2,  1;
       -2, -4, -2; 
        1,  2,  1]/4;

f22 =  [1, -2, 1;
        2, -4, 2; 
        1, -2, 1]/4;


f1m = reshape(f1(end:-1:1), [3,3]);
f2m = reshape(f2(end:-1:1), [3,3]);
f12m = reshape(f12(end:-1:1), [3,3]);
f11m = reshape(f11(end:-1:1), [3,3]);
f22m = reshape(f22(end:-1:1), [3,3]);

