import sncosmo
import triangle

class LC(object):
    """class to streamline light curve fits with SNCosmo """
    def __init__(self, model, data_file='sn.dat', truths=None):
        # super(, self).__init__()  <-- don't need this yet

        self.model = model
        self.data = sncosmo.read_lc(data_file)

        self._chi2 = None
        self._mcmc = None
        self._nest = None

    @property
    def chi2(self):
        return self.chi2

    @chi2.getter
    def chi2(self):
        if self._chi2 is None:
            print "running chi^2 fit"
            self._chi2 = sncosmo.fit_lc(self.data, self.model, vparam_names=['t0', 'x0', 'x1', 'c'],
                                        bounds={'c':(-0.3, 0.3), 'x1':(-3.0, 3.0)}, minsnr=3.0)
        return self._chi2

    @chi2.setter
    def chi2(self, val):
        self._chi2 = val
        return self._chi2

    @property
    def mcmc(self):
        return self.mcmc

    @mcmc.getter
    def mcmc(self):
        if self._mcmc is None:
            print "running MCMC fit"
            self._mcmc = sncosmo.mcmc_lc(self.data, self.model, vparam_names=['t0', 'x0', 'x1', 'c'],
                                        bounds={'c':(-0.3, 0.3), 'x1':(-3.0, 3.0)}, minsnr=3.0)
        return self._mcmc

    @mcmc.setter
    def mcmc(self, val):
        self._mcmc = val
        return self._mcmc

    @property
    def nest(self):
        return self.nest

    @nest.getter
    def nest(self):
        if self._nest is None:
            print "running nest fit"
            self._nest = sncosmo.nest_lc(self.data, self.model, vparam_names=['t0', 'x0', 'x1', 'c'],
                                    bounds={'c':(-0.3, 0.3), 'x1':(-3.0, 3.0)}, guess_amplitude_bound=True,
                                     minsnr=3.0, verbose=True)
        return self._nest

    @nest.setter
    def nest(self, val):
        self._nest = val
        return self._nest

#methods
#--------------

    def metadata(self):
        print self.model

    @staticmethod
    def statistics(self):
        print "chi2"
        print "dof: ", self._chi2[0].dof

    @staticmethod
    def plotLC(LC, fits=True):
        data = LC.data
        model = self.model
        chi2_model = self._chi2
        mcmc_model = self._mcmc
        nest_model = self._nest

        if fits:
            models=[model, chi2_model, mcmc_model, nest_model]
        else:
            models=[model]

        fig = sncosmo.plot_lc(data, models)

        return fig

    @staticmethod
    def plotCorner(self):
        mcmc = self._mcmc[0].vparam_names
        nest = self._nest[0].vparam_names

        mcmc_samples = self._mcmc[0].samples
        nest_samples = self._nest[0].samples

        mcmc_ndim, mcmc_nsamples = len(mcmc), len(mcmc.samples)

        nest_ndim, nest_nsamples = len(nest), len(nest.samples)

        #make figure

        figure_mcmc = triangle.corner(mcmc_samples, labels=[mcmc[0], mcmc[1], mcmc[2], mcmc[3]],
                         truths=[model.get(mcmc[0]), model.get(mcmc[1]),
                                 model.get(mcmc[2]), model.get(mcmc[3])],
                         range=mcmc_ndim*[0.9999],
                         show_titles=True, title_args={"fontsize": 12})

        figure_mcmc.gca().annotate("mcmc sampling", xy=(0.5, 1.0), xycoords="figure fraction",
                      xytext=(0, -5), textcoords="offset points",
                      ha="center", va="top")

        figure_nest = triangle.corner(nest_samples, labels=[nest[0], nest[1], nest[2], nest[3]],
                         truths=[model.get(nest[0]), model.get(nest[1]),
                                 model.get(nest[2]), model.get(nest[3])],
                         weights=self._nest[0].weights, range=nest_ndim*[0.9999],
                         show_titles=True, title_args={"fontsize": 12})

        figure_nest.gca().annotate("nest sampling", xy=(0.5, 1.0), xycoords="figure fraction",
                      xytext=(0, -5), textcoords="offset points",
                      ha="center", va="top")

        return figure_mcmc, figure_nest

    @staticmethod
    def plotTrace(self, mcmc=None, nest=None):
        if mcmc is None:
            mcmc = self._mcmc[0].vparam_names
            mcmc_samples = self._mcmc[0].samples
            mcmc_ndim, mcmc_nsamples = len(mcmc), len(mcmc_samples)

        else:
            mcmc = mcmc[0].vparam_names
            mcmc_samples = mcmc[0].samples
            mcmc_ndim, mcmc_nsamples = len(mcmc), len(mcmc_samples)

        if nest is None:
            nest = self._nest[0].vparam_names
            nest_samples = self._nest[0].samples
            nest_ndim, nest_nsamples = len(nest), len(nest_samples)


        trace_fig = plt.figure(figsize=(20,8))

        mcmc1 = trace_fig.add_subplot(241)
        mcmc2 = trace_fig.add_subplot(242)
        mcmc3 = trace_fig.add_subplot(243)
        mcmc4 = trace_fig.add_subplot(244)

        nest1 = trace_fig.add_subplot(245)
        nest2 = trace_fig.add_subplot(246)
        nest3 = trace_fig.add_subplot(247)
        nest4 = trace_fig.add_subplot(248)

        mcmc1.plot(mcmc_samples[:,0])
        mcmc2.plot(mcmc_samples[:,1])
        mcmc3.plot(mcmc_samples[:,2])
        mcmc4.plot(mcmc_samples[:,3])

        mcmc1.set_title('mcmc: ' + mcmc[0])
        mcmc2.set_title('mcmc: ' + mcmc[1])
        mcmc3.set_title('mcmc: ' + mcmc[2])
        mcmc4.set_title('mcmc: ' + mcmc[3])

        nest1.plot(nest_samples[:,0])
        nest2.plot(nest_samples[:,1])
        nest3.plot(nest_samples[:,2])
        nest4.plot(nest_samples[:,3])

        nest1.set_title('nest: ' + nest[0])
        nest2.set_title('nest: ' + nest[1])
        nest3.set_title('nest: ' + nest[2])
        nest4.set_title('nest: ' + nest[3])

        trace_fig.tight_layout()

        return trace_fig

    def comparefits2truth(self):
        pass

    def comparefits2fits(self):
        pass
