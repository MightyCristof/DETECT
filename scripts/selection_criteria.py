# scripts/selection_criteria.py
import pdb


#########################################
#
# Specific AGN Selection Criteria
#
#########################################
# Selection criteria must be a list of tuples (selection, operator, cutoff)
# Selection is either a color or magnitude (e.g. 'u-g', 'g-r', 'r', 'i')
# Operator is any in the set {'>=', '>', '<', '<='}
# Cutoff is the value of the 

analytic_criteria = [
    ('u-g','<',0.9),
    ('g-r','<',0.4),
]

richards2002_criteria = []


#########################################
#
# Main AGN Selection Class
#
#########################################
class AGNSelection:
    def __init__(self, criteria='analytic'):
        self.criteria = f"{criteria}_criteria"
        self.cuts = self._get_criteria()

    def _get_criteria(self):
        # Check if selection criteria exists in the current module
        global_vars = globals()
        if self.criteria in global_vars:
            return global_vars[self.criteria]
        else:
            raise ValueError(f"Unknown selection criteria: {self.criteria_name}")
    
    def is_agn(self, mags):
        selected = True
        for selection, operator, cutoff in self.cuts:
            if '-' in selection:
                color1, color2 = selection.split('-')
                svalue = mags[color1] - mags[color2]
            else:
                svalue = mags[selection]
            expression = f"{svalue} {operator} {cutoff}"
            # If any selection criteria are not met, return False
            if not eval(expression):
                selected = False
                break
         # If all selection criteria are met, return True
        return selected
