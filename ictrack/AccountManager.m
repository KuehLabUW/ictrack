classdef AccountManager
% Instantiate a listener for the bankaccount class
   methods (Static)      
      function assignStatus(BA)
         if BA.AccountBalance < 0  
            if BA.AccountBalance < -200
               BA.AccountStatus = 'closed';
            else
               BA.AccountStatus = 'overdrawn';
            end
         end
      end
      function addAccount(BA)
      % Add a listener for Insufficient Funds
         addlistener(BA, 'InsufficientFunds', ...
            @(src, evnt)AccountManager.assignStatus(src));
      end
   end
end
