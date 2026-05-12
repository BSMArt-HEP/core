"""

Helper functions for the xBit MLS scan


"""

import math
import numpy as np
import os
import random
from bsmart import core
import torch
import torch.nn as nn
import torch.optim as optim
from torch.autograd import Variable

"""

Functions from xBit aux

"""

def LINEAR_DIFF(a, b):
    return 1. * abs(a - b)

# gauss likelihood: do we need something else?
def gauss(val, m, s):
    return math.exp(-(val-m)**2/(2.*s)**2)


def logL(val, m, s):
    return 1./(1.+(val-m)**2/(2.*s)**2)


def limit_lh(steps):
    if (steps < 1.0e4):
        return 0.1
    elif (steps < 1.0e5):
        return 0.01
    elif (steps < 1.0e6):
        return 0.001
    else:
        return 1.0e-5


def step(x):
    if x > 0.5:
        return 1.
    else:
        return 0.


def exp_safe(x):
    try:
        out = math.exp(x)
    except OverflowError:
        out = float('inf')
    return out


def make_safe_lh(x):
    """
        curb extremes of a function, plateau at 100000
    """
    #return 100.0*math.tanh(x/100.0)
    return 1000000.0*math.tanh(x/1000000.0)

class NNMLS:
    def __init__(self, inputs):
        self.input_dim = len(inputs['Variables'])

        #self.train_lh = self.inputs['ML']['TrainLH']
        self.train_lh = eval(inputs['Setup'].get('TrainLH',"True"))


        
        if self.train_lh:
            self.log.info('Will train predictor on Likelihood')
            self.output_dim = 1
        else:
            self.log.info('Will train predictor on observables')
            #self.output_dim = len(self.inputs['Observables']) # this picks up all observables defined even in the tools
            # need to find how many observables we are *actually* using. 
            #self.output_dim = len(self.likelihoodfns)
            self.output_dim= len([x for x in self.observable_masks if x]) # these are false if masked
        
        self.log.info('Setting up predictor with %d outputs' % self.output_dim)
        self.epochs = inputs['Setup'].get('Epochs',500)
        self.iterations = inputs['Setup'].get('Iterations',20)
        self.neurons = inputs['Setup'].get('Neurons',[50,50,50])
        self.LR = inputs['Setup'].get('LR',0.001)
        
        if 'DensityPenalty' not in inputs['Setup']:
            self.distance_penalty = 0.1
        else:
            self.distance_penalty = inputs['Setup'].get('DensityPenality')  ## string because of eval
        
        
        self.use_classifier = eval(inputs['Setup'].get('Classifier',"False"))

    def get_train_lh(self):
        return self.train_lh

    def set_predictor(self):
        self.log.debug("Neurons: %s" %str(self.neurons))
        self.map1 = nn.Linear(self.input_dim, self.neurons[0])
        layers = [self.map1, nn.ReLU(), nn.Dropout(0.1)]

        for i in range(1,len(self.neurons)):
            layers.append(nn.Linear(self.neurons[i - 1], self.neurons[i]))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(0.1))

        layers.append(nn.Linear(self.neurons[-1], self.output_dim))

        self.predictor = nn.Sequential(*layers)
        self.predictor_optimizer = optim.Adam(self.predictor.parameters(),lr=self.LR)
        self.predictor_criterion = nn.MSELoss()

    def set_classifier(self):
        self.classifier = nn.Sequential(nn.Linear(self.input_dim, 50), nn.ReLU(), nn.Dropout(0.1),nn.Linear(50, 50), nn.ReLU(), nn.Dropout(0.1), nn.Linear(50,2), nn.Sigmoid())
        self.classifier_optimizer = optim.Adam(self.classifier.parameters(),lr=self.LR)
        self.classifier_criterion = nn.BCELoss()

    def splitter(self, x, y):   # split torch tensors into training and testing sets
        setlen=x.size(0)
        self.log.debug("Splitting dataset of size %d "%setlen)
        try:
            assert(y.size(0) == setlen)
        except Exception as e:
            self.log.info("Length of inputs/ouptuts %d vs  %d are not equal, %s" %(setlen,y.size(0),str(e)))
            raise
        perm = torch.randperm(setlen)
        splitind=int(0.8*float(setlen))
        idtrain = perm[:splitind]
        idtest=perm[splitind:]
        return x[idtrain], y[idtrain], x[idtest], y[idtest]

    def do_train(self, x_val, y_val, name):
        #self.log.info(f"calling do_train with {len(x_val)} points and {self.inputs['Setup']['Cores']} cores")
        self.log.info("calling do_train with %d points and %d cores" %(len(x_val),self.inputs['Setup']['Cores']))
        #print(x_val)
        #print(y_val)
        #print(f"x: {x_val} \n y: {y_val}")
        #self.classifier = nn.Sequential(nn.Linear(self.input_dim, 50), nn.ReLU(), nn.Dropout(0.1),nn.Linear(50, 50), nn.ReLU(), nn.Dropout(0.1), nn.Linear(50,2), nn.Sigmoid())
        torch.set_num_threads(self.inputs['Setup']['Cores'])
        
        best = 1.0e6        
        wait = 0
        batch_size = 64
        nr_batches = max(1, int(len(x_val) / batch_size))
        #patience = max(100, self.epochs / 5 / nr_batches)
        patience = max(20, self.epochs / 5 / nr_batches)  ## number of epochs to wait if there is no improvement in loss

        #self.log.info(f"Patience is {patience}")
        self.log.info("Patience is %.1f" %patience)
        #torchx=torch.tensor(x_val)
        #torchx=x_val

        if name is self.predictor:
            criterion = nn.MSELoss()                                        #self.predictor_criterion
            optimizer = optim.Adam(self.predictor.parameters(),lr=self.LR)  #self.predictor_optimizer
            torchx=x_val # feed it a torch tensor if it is a predictor
            torchy=torch.tensor(y_val).reshape((-1,self.output_dim))
            checkpoint_name=os.path.join(self.inputs['Setup']['RunName'],'checkpoint_predictor.pt')
        else:
            criterion = nn.BCELoss()                                        #self.classifier_criterion
            optimizer = optim.Adam(self.classifier.parameters(),lr=self.LR) #self.classifier_optimizer
            torchx=torch.tensor(x_val)
            torchy=torch.tensor(y_val).reshape((-1,2))  # for some stupid reason have two outputs
            checkpoint_name=os.path.join(self.inputs['Setup']['RunName'],'checkpoint_classifier.pt')
        
        
        x_train, y_train, x_test, y_test = self.splitter(torchx,torchy)

        permutation = torch.randperm(x_train.size(0))

        ## want to get the best lost from *previous* run!
        if os.path.exists(checkpoint_name):
            old_checkpoint = torch.load(checkpoint_name)
            previous_best = old_checkpoint['loss']

        for epoch in range(self.epochs):
            name.train()
            for _ in range(x_train.size()[1]):
                for i in range(0, x_train.size(0), batch_size):
                    optimizer.zero_grad()
                    indices = permutation[i:i+batch_size]
                    x_batch = x_train[indices]
                    y_batch = y_train[indices] #[indices]
                    self.loss = criterion(name(Variable(x_batch)), Variable(y_batch))
                    self.loss.backward()
                    optimizer.step()
                name.eval()

            y_pred = name(Variable(x_test).data.cpu())
            loss = criterion(Variable(y_pred), Variable(y_test))

            if epoch == 0:  # have to deal with the case that the initial loss is *worse* than 1e6!
                best = loss

            if epoch % int(self.epochs/10) == 0:
                self.log.info("Epoch: %i; loss %f, best loss %f" %(epoch, loss,best))

            if loss < best:
                best = loss
                wait = 0
                
                torch.save({'state_dict': name.state_dict(), 'optimizer': optimizer.state_dict(), 'loss': self.loss}, checkpoint_name)
            else:
                wait = wait + 1

            if (wait > patience or epoch == (self.epochs - 1)) and os.path.exists(checkpoint_name):
                """
                    This stops training when there is no more gain and reloads the best network
                """
                self.log.info("Loading state dict %s" % checkpoint_name)
                checkpoint = torch.load(checkpoint_name)
                name.load_state_dict(checkpoint['state_dict'])
                optimizer.load_state_dict(checkpoint['optimizer'])
                self.loss = checkpoint['loss']
                break
        #self.log.info(f"Completed training after {epoch+1} epochs")
        self.log.info("Completed training after %d epochs" % (epoch+1))

    def scale_data(self, input):
        """

            This is needed to curb excess values 
            'input' is a list of observables
       
        """
        self.scalings = []
        transposedinput = np.array(input).T
        data_out = None
        for data in transposedinput:
            #print(data)
            if max(data) / min(data) > 100:
                if data_out is None:
                    data_out = np.array([np.log(data)])
                else:
                    data_out = np.concatenate((data_out, [np.log(data)]), axis = 0)
                self.scalings.append("log")
            else:
                if data_out is None:
                    data_out = np.array([data])
                else:
                    data_out = np.concatenate((data_out, [data]), axis = 0)
                self.scalings.append("id")

        return((data_out.T).tolist())
    
    def rescale_data(self, input):
        """
            undoes the effect of the rescaling, to return to the actual observables
        """
        transposedinput = np.array(input).T
        data_out = None
        for data, scale in zip(transposedinput, self.scalings):
            if scale == "log":
                if data_out is None:  # why doesn't it work with empty lists?
                    data_out = np.array([np.exp(data)])
                else:
                    data_out = np.concatenate((data_out, [np.exp(data)]),axis=0)
            elif scale == "id":
                if data_out is None:
                    data_out = np.array([(data)])
                else:
                    data_out = np.concatenate((data_out, [data]), axis=0)
        return((data_out.T).tolist())