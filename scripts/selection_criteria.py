# scripts/selection_criteria.py
import pdb


#########################################
#
# Main AGN Selection Class
#
#########################################
class AGNSelection:
    def __init__(self, criteria='Analytic'):
        self.criteria = f"{criteria}Criteria"
        self.func = self._get_criteria()

    def _get_criteria(self):
        # Check if selection criteria exists in the current module
        global_vars = globals()
        if self.criteria in global_vars:
            criteria_class = global_vars[self.criteria]
            if callable(criteria_class):
                return criteria_class()
            else:
                raise TypeError(f"{self.criteria_name} is not callable.")
        else:
            raise ValueError(f"Unknown selection criteria: {self.criteria_name}")

    def is_agn(self, mags):
        return self.func.select_agn(mags)


#########################################
#
# Specific AGN Selection Classes
#
#########################################
class AnalyticCriteria:
    def __init__(self):
        self.cuts = {
            'u-g': 0.9,
            'g-r': 0.4
        }

    def select_agn(self, mags):
        if (mags['u']-mags['g'] < self.cuts['u-g']) and (mags['g']-mags['r'] < self.cuts['g-r']):
            return True
        else:
            return False


class Richards2002Criteria:
    def __init__(self):
        self.cuts = {}

    def select_agn(self):
        # Implement selection for Richards+2002
        print("AGN selected using Richards+2002 criteria")
