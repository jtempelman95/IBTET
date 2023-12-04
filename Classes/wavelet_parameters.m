classdef wavelet_parameters

    properties
        fo % Center freq
        fl % lower freq
        fu % upper freq
        Fs % sampling rate
        fn %N Freq bins
    end

    methods
        function self = wavelet_parameters(fo,fl,fu,Fs,fn)
         if nargin>0
            self.fo = fo;
         end
         if nargin > 1
            self.fl = fl;
         end
         if nargin > 2
            self.fu = fu;
         end
         if nargin >3
            self.Fs = Fs;
         end
         if nargin >4
            self.fn = fn;
         end
         end
    end

end