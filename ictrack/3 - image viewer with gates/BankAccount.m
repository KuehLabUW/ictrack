classdef BankAccount < handle
   properties (Hidden)
      AccountStatus = 'open';
   end
   
   properties (SetAccess = private)
      AccountNumber
      AccountBalance = 0;
   end
   
   events
      InsufficientFunds
   end
   
   methods
      function BA = BankAccount(AccountNumber,InitialBalance)
         if nargin < 2
            error('BankAccount:InvalidInitialization',...
               'Must provide an account number and initial balance')
         end
         BA.AccountNumber  = AccountNumber;
         BA.AccountBalance = InitialBalance;
         %AccountManager.addAccount(BA);
      end % BankAccount
      
      function deposit(BA,amt)
         BA.AccountBalance = BA.AccountBalance + amt;
         if BA.AccountBalance > 0
            BA.AccountStatus = 'open';
         end
      end % deposit
      
      function withdraw(BA,amt)
         if (strcmp(BA.AccountStatus,'closed')&& BA.AccountBalance < 0)
            disp(['Account ',num2str(BA.AccountNumber),' has been closed.'])
            return
         end
         newbal = BA.AccountBalance - amt;
         BA.AccountBalance = newbal;
         if newbal < 0
            notify(BA,'InsufficientFunds')
         end
      end % withdraw
      
   end % methods
end % classdef