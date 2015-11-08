import sncosmo
import triangle
import matplotlib.pyplot as plt

class LC(object):
    """class to streamline light curve fits with SNCosmo """
    def __init__(self, model, data_file='sn.dat', truths=None):
        # super(, self).__init__()  <-- don't need this yet

        self.model = model
        self.data = sncosmo.read_lc(data_file)

        # no host extinction
        self._fitOut = None
        self.fitRes = None
        self.fitModel = None

        self._mcmcOut = None
        self.mcmcRes = None
        self.mcmcModel = None

        self._nestOut = None
        self.nestRes = None
        self.nestModel = None

        # with host extinction
        self._fitExtOut = None
        self.fitExtRes = None
        self.fitExtModel = None

        self._mcmcExtOut = None
        self.mcmcExtRes = None
        self.mcmcExtModel = None

        self._nestExtOut = None
        self.nestExtRes = None
        self.nestExtModel = None


    @property
    def fitOut(self):
        return self.fitOut

    @fitOut.getter
    def fitOut(self):
        if self._fitOut is None:
            print "running chi^2 fit"

            vparams = ['t0', 'x0', 'x1', 'c']

            self._fitOut = self.runFit(vparams)
            self.fitRes = self._fitOut[0]
            self.fitModel = self._fitOut[1]

        return self._fitOut

    @fitOut.setter
    def fitOut(self, val):
        self._fitOut = val
        self.fitRes = self._fitOut[0]
        self.fitModel = self._fitOut[1]

        return self._fitOut

    @property
    def mcmcOut(self):
        return self.mcmcOut

    @mcmcOut.getter
    def mcmcOut(self):
        if self._mcmcOut is None:
            print "running MCMC fit"

            vparams = ['t0', 'x0', 'x1', 'c']

            self._mcmcOut = self.runMCMC(vparams)
            self.mcmcRes = self._mcmcOut[0]
            self.mcmcModel = self._mcmcOut[1]

        return self._mcmcOut

    @mcmcOut.setter
    def mcmcOut(self, val):
        self._mcmcOut = val
        self.mcmcRes = self._mcmcOut[0]
        self.mcmcModel = self._mcmcOut[1]

        return self._mcmcOut

    @property
    def nestOut(self):
        return self.nestOut

    @nestOut.getter
    def nestOut(self):
        if self._nestOut is None:
            print "running nest fit"

            vparams = ['t0', 'x0', 'x1', 'c']

            self._nestOut = self.runNest(vparams)
            self.nestRes = self._nestOut[0]
            self.nestModel = self._nestOut[1]

        return self._nestOut

    @nestOut.setter
    def nestOut(self, val):
        self._nestOut = val
        self.nestRes = self._nestOut[0]
        self.nestModel = self._nestOut[1]

        return self._nest


#methods
#--------------

    def runFit(self, vparams):
        fitOut = sncosmo.fit_lc(self.data, self.model, vparam_names=vparams,
                                    bounds={'c':(-0.3, 0.3), 'x1':(-3.0, 3.0)}, minsnr=3.0)
        return fitOut

    def runMCMC(self, vparams):
        mcmcOut = sncosmo.mcmc_lc(self.data, self.model, vparam_names=vparams,
                                    bounds={'c':(-0.3, 0.3), 'x1':(-3.0, 3.0)}, minsnr=3.0)
        return mcmcOut

    def runNest(self, vparams):
        nestOut = sncosmo.nest_lc(self.data, self.model, vparam_names=['t0', 'x0', 'x1', 'c'],
                                bounds={'c':(-0.3, 0.3), 'x1':(-3.0, 3.0)}, guess_amplitude_bound=True,
                                 minsnr=3.0, verbose=True)
        return nestOut

    def metadata(self):
        print self.model

    @staticmethod
    def statistics(self):
        print "chi2"
        print "dof: ", self._fitOut[0].dof

    @staticmethod
    def plotLC(LC , fits=True):
        data = LC.data
        model = LC.model
        fit_model = LC.fitModel
        mcmc_model = LC.mcmcModel
        nest_model = LC.nestModel

        if fits:
            models=[model, fit_model, mcmc_model, nest_model]
        else:
            models=[model]

        fig = sncosmo.plot_lc(data, model=models)

        return fig

    @staticmethod
    def plotCorner(LC):
        model = LC.model

        mcmcVParams = LC.mcmcRes.vparam_names
        nestVParams = LC.nestRes.vparam_names

        mcmcSamples = LC.mcmcRes.samples
        nestSamples = LC.nestRes.samples

        mcmc_ndim, mcmc_nsamples = len(mcmcVParams), len(mcmcSamples)
        nest_ndim, nest_nsamples = len(nestVParams), len(nestSamples)

        #make figure

        figure_mcmc = triangle.corner(mcmcSamples, labels=[mcmcVParams[0], mcmcVParams[1], mcmcVParams[2], mcmcVParams[3]],
                         truths=[model.get(mcmcVParams[0]), model.get(mcmcVParams[1]),
                                 model.get(mcmcVParams[2]), model.get(mcmcVParams[3])],
                         range=mcmc_ndim*[0.9999],
                         show_titles=True, title_args={"fontsize": 12})

        figure_mcmc.gca().annotate("mcmc sampling", xy=(0.5, 1.0), xycoords="figure fraction",
                      xytext=(0, -5), textcoords="offset points",
                      ha="center", va="top")

        figure_nest = triangle.corner(nestSamples, labels=[nestVParams[0], nestVParams[1], nestVParams[2], nestVParams[3]],
                         truths=[model.get(nestVParams[0]), model.get(nestVParams[1]),
                                 model.get(nestVParams[2]), model.get(nestVParams[3])],
                         weights=LC.nestRes.weights, range=nest_ndim*[0.9999],
                         show_titles=True, title_args={"fontsize": 12})

        figure_nest.gca().annotate("nest sampling", xy=(0.5, 1.0), xycoords="figure fraction",
                      xytext=(0, -5), textcoords="offset points",
                      ha="center", va="top")

        return figure_mcmc, figure_nest

    @staticmethod
    def plotTrace(LC):
        mcmcVParams = LC.mcmcRes.vparam_names
        nestVParams = LC.nestRes.vparam_names

        mcmcSamples = LC.mcmcRes.samples
        nestSamples = LC.nestRes.samples

        trace_fig = plt.figure(figsize=(20,8))

        mcmc1 = trace_fig.add_subplot(241)
        mcmc2 = trace_fig.add_subplot(242)
        mcmc3 = trace_fig.add_subplot(243)
        mcmc4 = trace_fig.add_subplot(244)

        nest1 = trace_fig.add_subplot(245)
        nest2 = trace_fig.add_subplot(246)
        nest3 = trace_fig.add_subplot(247)
        nest4 = trace_fig.add_subplot(248)

        mcmc1.plot(mcmcSamples[:,0])
        mcmc2.plot(mcmcSamples[:,1])
        mcmc3.plot(mcmcSamples[:,2])
        mcmc4.plot(mcmcSamples[:,3])

        mcmc1.set_title('mcmc: ' + mcmcVParams[0])
        mcmc2.set_title('mcmc: ' + mcmcVParams[1])
        mcmc3.set_title('mcmc: ' + mcmcVParams[2])
        mcmc4.set_title('mcmc: ' + mcmcVParams[3])

        nest1.plot(nestSamples[:,0])
        nest2.plot(nestSamples[:,1])
        nest3.plot(nestSamples[:,2])
        nest4.plot(nestSamples[:,3])

        nest1.set_title('nest: ' + nestVParams[0])
        nest2.set_title('nest: ' + nestVParams[1])
        nest3.set_title('nest: ' + nestVParams[2])
        nest4.set_title('nest: ' + nestVParams[3])

        trace_fig.tight_layout()

        return trace_fig

    def compareExtinction(self):
        figCompExt = plt.figure()

        fitComp = figCompExt.add_subplot(131)
        mcmcComp = figCompExt.add_subplot(132)
        nestComp = figCompExt.add_subplot(133)

        fitComp.plot_lc()

    def comparefits2truth(self):
        pass

    def comparefits2fits(self):
        pass

    def calculateBias(LC):
        model = LC.model
        mcmcSamples = LC.mcmcRes.samples
        nestSamples = LC.nestRes.samples

        #fitBias = [mcmcSamples[0,x] - model.get() for x in ]


    def calculateVariance(self):
        pass

    def compFitExt():
        pass

    def compMcmcExt():
        pass

    def compNestExt():
        pass
