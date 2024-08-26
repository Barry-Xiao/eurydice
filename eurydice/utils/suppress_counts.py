from functools import partial
import itertools as it

import numpy as np
import pandas as pd


def suppress_counts(counts, low_thresh=5, axis=-1, percent_axis=None, 
                    total=False, round_up=False, squish_cols=None, 
                    penalty_rows=None, 
                    optimize='min_counts', sum_min_scale=50):
    """
    Suppress counts below the threshhold for display

    ... 

    Parameters
    ----------
    counts: pd.DataFrame
        The counts to be suppressed. We semi assume that cohorts are the
        rows and groupings are the columns.
    low_thresh: float
        Any counts below this value should be suppressed.
    axis: {0, 1, -1, 'index', 'columns', 'both'}
        The axis over which to sum and suppress
    percent_axis: {0, 1, -1, 'index', 'columns', 'both', None}
        The axis to use for calculating the percent. If not specified, this
        defaults to the axis used for suppression.
    total: bool, optional
        Should the row and column totals be shown
    round_up: bool, optional
        If a suppressed value is a multiple fo the low threshhold, should the
        masked value be rounded up for more accurate suppression
    squish_cols: str, list, optional
        A column or list of columns in counts that should be used to
        absorb supression first, where possible. In two dimesnional
        data, these will be columns. The columsn will be preferentially
        ordered based on what is presented in the data
    penalty_rows: str, list
        The rows which should not be used to absorb values, unless
        necessary necessary. If the data is two dimensional, these will
        be rows. For example, if a cohort has a known sum and only one
        group, it may be undesirable to suppress those values.
    optimize: {'min_counts', 'overlap_min'}
        How the secondary suppression should be optimized, if there are
        multiple candiate values.
            - "min_counts" will suppress the smallest counts in a given
                row or column, even if that means suppressing multiple
                values same row or column.
            - "overlap_min" will try to minimize the number of cells per
                row and the minimum sum by adding a count penalty (set
                with `min_sum_scale` to the cells with multiple values.
    min_sum_scale: int, optional
        A parameter to help tune the balance between minimizing the cells
        per row and minimizing the suppressed counts in those cells.
        Used only for "overlap_min" optimization.

    Returns
    -------
    pd.DataFrame
        The suppressed counts with the masked values rounded to allow
        suppression
    pd.DataFrame
        The boolean mask that shows which cells have been suppressed
    pd.DataFrame
        The suppressed percent with the masked values rounded for
        suppression

    Raises
    ------
    ValueError
        The axis value is misspecified
    ValueError
        If the optimization isnt solved correctly
    """

    # Standardizes the axis arguement. If columns is specified, we just
    # transpose the data and then operate on the row data.
    if axis in {0, '0', 'columns'}:
        axis = 0
        counts = counts.T
    elif axis in {1, 'index', '1'}:
        axis = 1
    elif axis in {-1, 'both', '-1', None}:
        axis = -1
    else:
        raise ValueError('The axis value is misspectified. Possible '
                         'values are "index" (1), "columns" (0), and '
                         '"both" (-1)')
    if percent_axis is None:
        percent_axis = axis
        
    # ########################### Total values ###########################
    # Calculates row and sum values
    rows = counts.sum(axis=1)
    cols = counts.sum(axis=0)
    totals = rows.sum()

    # Finds rows and columns that need to be suppressed
    rows_below = (rows < low_thresh) & (rows > 0)
    cols_below = (cols < low_thresh) & (cols > 0)
        
    # ########################### General masking ########################
    # Masks the vlaues based on the count values, identifying the 
    # positions that are zeros, below the threshhold, and ≥threshhold
    # The primary suppression is relatively easy: we're looking at values
    # below the threshhold
    zeros = (counts == 0)
    below = (counts < low_thresh) & ~(zeros)
    above = ~(zeros | below)

    # Reference matrix 
    blank = counts.isna()

    # ####################################################################
    # ##################### Secondary Suppression Masks ##################
    # ####################################################################

    # ##################### Specialty suppression masks ##################
    # If we're interested in using a specific column to absorb 
    # suppression, we want to be able to designate that
    squish_mask, squish_cols, squish_order = \
        _build_special_suppression(squish_cols, axis, blank)
    squish_mask = squish_mask & above

    # We only penalize rows currently
    penalty_mask, penalty_rows, _ = \
        _build_special_suppression(penalty_rows, 1, blank)

    # We want a suppression mask for the row total and column total, in 
    # case we need to suppress the sums
    if axis < 0:
        row_suppress, row_sum_below = _build_sum_suppression(
            rows=rows, 
            low_thresh=low_thresh,
            squish_row=blank.any(axis=1),
            penalty_row=penalty_mask.any(axis=1),
            squish_order=blank.any(axis=1) * 0 + 99,
            row_count=(~zeros).sum(axis=1),
            )
        col_suppress, col_sum_below = _build_sum_suppression(
            rows=cols,
            low_thresh=low_thresh,
            squish_row=squish_mask.any(axis=0),
            penalty_row=blank.any(axis=0),
            squish_order=squish_order,
            row_count=(~zeros).sum(axis=0),
            )
    else:
        row_suppress, row_sum_below = _build_sum_suppression(
            rows=rows,
            low_thresh=low_thresh,
            squish_row=squish_mask.any(axis=1),
            penalty_row=penalty_mask.any(axis=1),
            squish_order=squish_order,
            row_count=(~zeros).sum(axis=1),
            )
        col_suppress = blank.all(axis=0)
        col_sum_below = blank.all(axis=0)

    # Secondary suppression is an optimization problem where we're trying
    # to pick the best cells to suppress so that every cell that needs to
    # be suppressed is suppressed but we dont suppress any more cells than
    # need to be suppressed. 
    # 
    # We do this by finding every place that needs to be suppressed and
    # optimizing the rows and columns that can be used to mask the data.
    # Then we check that the rows and columns meet out suppression criteria
    # (basically the optimal mask has to actually suppress everything).
    # If that doens't work, we incldue the cell which is unoptimized in
    # as a below cell and then we try re-optimizing. 
    #
    # As a programming note, I know "where" loops are inelegant, but I
    # dont actually think there's a terribly elegant way to do this. 
    row_optimum = blank.copy()
    col_optimum = blank.copy()

    solution, suppress_match_r, suppress_match_c = \
        _check_overall(below, zeros, counts, axis, low_thresh)

    for i in range(5):
        # We find the missing cells and add those to the below matrix
        # and just treat them as below the threshhold
        col_missed = \
            (~suppress_match_r.to_frame().values & col_optimum)
        row_missed = \
            (~suppress_match_c.to_frame().values.T & row_optimum)

        # If there are multiple canidate values in the same column, pick the
        # smallest value
        if (col_missed.sum(axis=0) > 1).any():
            cmiss_vals = counts.copy().mask(~col_missed)
            col_missed = (cmiss_vals == cmiss_vals.min(axis=0).values.T) * \
                col_missed
        # Multiple canidate values in a row? Pick the minimum. 
        if (row_missed.sum(axis=1) > 1).any():
            rmiss_val = counts.copy().mask(~row_missed)
            row_missed = \
                (rmiss_val == rmiss_val.min(axis=1).to_frame().values) & \
                row_missed

        below = below | row_missed | col_missed
        above = above & ~below


        row_below, col_below, row_nz, col_nz, row_sum_supp, \
                col_sum_supp, canidates =  \
            _calculate_mask_properties(below, above, blank, zeros, 
                                        row_suppress, col_suppress)

        # And then check for new optimized data
        row_maximize, row_optimum = _optimize_row_suppression_mask(
            row_below=row_below,
            col_below=col_below,
            squish_mask=squish_mask & above,
            squish_order=squish_order,
            row_sum_supp=row_sum_supp,
            col_sum_supp=col_sum_supp,
            col_nz=col_nz,
            penalty_mask=penalty_mask & (axis >= 0) & above,
            canidates=canidates,
            axis=axis,
            )
        col_maximize, col_optimum = _optimize_row_suppression_mask(
            row_below=col_below.T,
            col_below=row_below.T & (axis < 0),
            squish_mask=blank.T,
            squish_order=blank.any(axis=1) * 1 + 98,
            row_sum_supp=col_sum_supp.T,
            col_sum_supp=row_sum_supp.T,
            col_nz=row_nz.T,
            penalty_mask=penalty_mask.T,
            canidates=canidates.T,
            axis=-1,
            )
        col_optimum = col_optimum.T & (axis < 0)
        col_maximize = col_maximize.T

        # Checks if the mask solves the issue; if we can solve, we continue
        solution, suppress_match_r, suppress_match_c = _check_overall(
            mask=(col_optimum | row_optimum | below),
            zeros=zeros,
            counts=counts,
            row_suppress=(row_suppress | row_sum_below) & total,
            col_suppress=(col_suppress | col_sum_below) & total,
            axis=axis,
            low_thresh=low_thresh,
            )

        if solution:
            break
    if ~solution:
        raise ValueError('The mask was unsolved')

    # print(row_optimum * 1)
    # print(col_optimum * 1)
    # print(below * 1)
        
    # ######################### Solve the Optimization ###################
    if (~row_optimum & ~col_optimum).all().all():
        # If there's nothing to optimize, we take the results
        masked = (below | row_sum_supp | col_sum_supp) & ~zeros
    elif (row_optimum.sum(axis=1) < 2).all() & \
            (col_optimum.sum(axis=0) < 2).all():
        masked = (below | row_optimum |  col_optimum) & ~zeros
    else:        
        # At this point, we're going to essentially take the optimum rows and
        # columns and build a set of masks that satisify our critiera. So, we
        # start by getting a listing of the cells that are below the
        # suppression limit
        below_long = _list_below(below)

        # For each cells, we construct a set of possible masks
        mask_list = [
            _make_coord_mask(row, col, below, row_optimum, col_optimum)
            for _, [row, col] in below_long[['index', 'columns']].iterrows()
            ]

        # And then we calculate the possible combinations of thoses masks and
        # propreties of the masks. In addition to existing masks that solve
        # the optimization, we'll also suppress any non-zero value in a row
        # or column that requires sum suppression.
        # We then calculate a set of parameters for the mask.
        sum_mask = (row_sum_supp | col_sum_supp) & ~zeros

        mask_check = partial(
            _check_mask_properties,
            counts=counts,
            sum_mask=sum_mask & total,
            zeros=zeros,
            below=below,
            axis=axis,
            sum_min_scale=sum_min_scale,
            low_thresh=low_thresh,
            row_suppress=(row_suppress | row_sum_below) & total,
            col_suppress=(col_suppress | col_sum_below) & total,
            )
        
        mask_summaries = pd.DataFrame([
            mask_check(masks=m) for m in it.product(*mask_list)
            ])
        mask_summaries = \
            mask_summaries.loc[mask_summaries['solution']].copy()
        mask_summaries.drop_duplicates(
            subset=['num_cells', 'min_counts', 'min_cells', 'min_overlap', 
                    'solution'],
            inplace=True
            )
        mask_summaries.sort_values(
            by=['num_cells', 'min_counts', 'min_cells', 'min_overlap'],
            ascending=True, 
            inplace=True)
        # And then we sort by the mode and pull out the best mask that solves
        # things
        mask_summaries.sort_values(['num_cells', optimize],
                                   ascending=True,
                                   inplace=True)
        masked = mask_summaries.iloc[0]['mask']
    
    # ####################################################################
    # ##################### Masking #####################
    # ####################################################################
    # And now we start masking. We'll first suppress the counts that 
    # need to be suppressed 
    round_check = \
        (np.ceil(counts / low_thresh) == (counts / low_thresh)) * \
            masked * round_up
    counts_supp = \
        counts.mask(masked, 
                    np.ceil((counts + round_check) / low_thresh) * low_thresh)

    # ### TOTALS ###
    # We add in the suppressed totals
    if total & (axis in {-1}):
        row_suppress = (row_suppress | row_sum_below)
        col_suppress = (col_suppress | col_sum_below)
        
        counts_supp['total'] = \
            rows.mask(row_suppress, np.ceil(rows / low_thresh) * low_thresh)
        counts_supp.loc['total', cols.index] = \
            cols.mask(col_suppress, np.ceil(cols / low_thresh) * low_thresh)
        counts_supp.loc['total', 'total'] = totals
        masked.loc['total', cols.index] = col_suppress
        masked['total'] = row_suppress
        masked.loc['total', 'total'] = False
    elif total:
        row_suppress = (row_suppress | row_sum_below)
        counts_supp['total'] = \
            rows.mask(row_suppress, np.ceil(rows / low_thresh) * low_thresh)
        masked['total'] = row_suppress

    masked = (masked == True)
        
    # ### Calculates the percent mask ###
    if percent_axis == -1:
        perc_supp = counts_supp / totals
    elif percent_axis == 0:
        perc_supp = counts_supp.div(cols, axis=1)
    else:
        perc_supp = counts_supp.div(rows, axis=0)

    perc_supp.mask(perc_supp > 1, 1, inplace=True)

    if axis == 0:
        counts_supp = counts_supp.T
        masked = masked.T
        perc_supp = perc_supp.T

    return counts_supp.astype(int), masked, perc_supp


def _build_special_suppression(cols, axis, blank):
    """
    Builds a specialty suppresion mask for squish columsn or penalty rows
    """
    # Gets an empty array (all false)
    mask_ = blank & False
    if cols is None:
        return mask_, [], mask_.any(axis=(axis >= 0) * 1) * 0 + 99

    # Converts the columns to a list
    if isinstance(cols, (str, float, int)):
        cols = [cols]

    # If we're looking at a 2D suppression, we want suppress columns. 
    # Otherwise, we suppress the rows
    if axis < 0:
        supp_cols = [c for c in cols if c in blank.columns]
        if len(supp_cols) > 0:
            mask_[supp_cols] = True
        order_ = pd.Series(data=np.arange(len(supp_cols))[::-1],
                           index=supp_cols)
        supp_order = (mask_.any(axis=0) + order_).fillna(99)        
    else:
        supp_cols = [c for c in cols if c in blank.index]
        if len(supp_cols) > 0:
            mask_.loc[supp_cols, mask_.columns] = True
        order_ = pd.Series(data=np.arange(len(supp_cols)),
                           index=supp_cols)
        supp_order = (mask_.any(axis=1) + order_).fillna(99)

    return mask_, supp_cols, supp_order


def _build_sum_suppression(rows, low_thresh, squish_row, penalty_row, squish_order, row_count, other_supp=False):
    """
    Identify row (or column) sums that need secondary suppression
    """
    # We find the number of row that are below the threshhold
    rows_below = ((rows < low_thresh) & (rows > 0)) | other_supp
    
    if rows_below.sum() != 1:
        # If there are no rows or more than 1 row below the threshhold, 
        # we just let htem suppress eachother and contineu with our lives
        row_suppress = rows_below & False
        return row_suppress, rows_below

    elif (squish_row & ~rows_below).any():
        # If there's one row below and an eligible sqiush row (i.e. the 
        # low row is not the squish row), than we take the first squish 
        # row that was specified

        squish_min = squish_order.mask(rows_below, 99)
        squish_min = squish_min == squish_min.min()
        
        return squish_min, rows_below
    else:
        # Otherwise, we start working through the remaining options. 
        
        # We start with the number of non-zero values in the row. I want 
        # to minimize the impact, so the goal is first to suppress groups
        # with only one cell filled; and then second to suppress columns
        # where he total counts are minimized.

        # So, we start by adding a count penalty to the penalty mask:
        # anything that should be penalized will be given a non-zero
        # count greater than the maximum non zero count currently in the
        # table, so it gets ignored
        row_count.mask(penalty_row, penalty_row.max() + 1, inplace=True)
        row_count.mask(rows_below, np.nan, inplace=True)

        if (row_count == 1).sum() == 1:
            # If there's only one row that has one value, we choose that
            # one
            return (row_count == 1), rows_below
        elif (row_count == 1).sum() > 1:
            # if we have more than one row which has one value, we find
            # the minimum row
            row_min = rows.copy().mask(rows_below)
            row_min = row_min.mask((row_count != 1)).min()
            return (rows == row_min), rows_below
        else:
            # Otherwise, we need to penalize the penalty mask and then
            # we'll take the minimum value 
            row_min = rows.copy().mask(penalty_row, rows + rows.max() / 2)
            row_min = row_min.mask(rows_below).min()
            return (rows == row_min), rows_below


def _check_mask_properties(masks, counts, sum_mask, zeros, below, row_suppress, col_suppress, axis, sum_min_scale=100, low_thresh=5):
    """
    Builds an overlapping mask and checks the mask properties
    """
    # Builds the overall mask
    combined = pd.concat(axis=1, objs=[m.unstack() for m in masks])
    combined = combined.any(axis=1).unstack().T
    mask = below | combined | sum_mask
    
    # Checks if the mask is a valid mask    
    is_solution, _, _ = _check_overall(mask, zeros, counts, axis, low_thresh)
    
    # If the mask is a solution, we want to calculate the mask sum, the number
    # of cells suppressed, and the penalized overlap sum
    if is_solution:
        # Calculates the min overlap
        min_overlap = (mask * counts) + \
            (mask * (mask.sum(axis=1).to_frame().values - 1)) * sum_min_scale + \
            (mask * (mask.sum(axis=0).to_frame().values.T - 1)) * sum_min_scale
        num_rows = (mask * 1).sum(axis=1)
        num_cols = (mask * 1).sum(axis=0)
        solution = {
            'num_cells': (mask * 1).sum().sum(),
            'min_cells': num_rows.mask(num_rows == 0).mean() + \
                num_cols.mask(num_cols == 0).mean(),
            'min_counts': (mask * counts).sum().sum(),
            'min_overlap': min_overlap.sum().sum(),
            'mask': mask,
            'solution': is_solution,
            }
    else:
        solution = {
                   'num_cells': np.nan,
                   'min_counts': np.nan,
                   'min_cells': np.nan,
                   'min_overlap': np.nan,
                   'mask': np.nan,
                   'solution': is_solution,
                   }
    return pd.Series(solution)


def _calculate_mask_properties(below, above, blank, zeros, row_suppress, col_suppress):
    """
    Calculates the row properties for secondary suppression
    """
    # Looks for values are canidates for suppression (above the threshold
    # and in a row or column with a value that needs to be suppressed)
    row_below = ((below * 1).sum(axis=1) == 1).to_frame().values & above
    col_below = ((below * 1).sum(axis=0) == 1).to_frame().values.T & above

    # Counts the number of non-zero values per row and column
    row_nz = \
        (((~zeros).sum(axis=1).to_frame().values + blank * 0) > 1) & above
    col_nz = \
        (((~zeros).sum(axis=0).to_frame().values.T + blank * 0) > 1) & above

    # Identifies the rows and columns that need to be fully suppressed?
    row_sum_supp = above & row_suppress.to_frame().values
    col_sum_supp = above & col_suppress.to_frame().values.T

    # Identifies cells avaliable for secondary masking
    canidates = \
        row_below | col_below | row_sum_supp | col_sum_supp & above       

    res = (row_below, col_below, row_nz, col_nz, row_sum_supp, 
           col_sum_supp, canidates)

    return res


def _check_overall(mask, zeros, counts, row_suppress, col_suppress, axis=-1, low_thresh=5):
    """
    Checks a mask suppresses everything that needs to be suppressed

    We use the critiera that for a given row or column, there are three
    valid options for suppression:
    1. No values are suppressed
    2. At least 2 values are suppressed
    3. If is one value suppressed, everything else is 
    """
    suppressed = ((counts.mask(mask) >= low_thresh) | zeros | mask)

    # Rows
    row_check = (mask.sum(axis=1) == 0) | (mask.sum(axis=1) >= 2) | \
        ((mask | zeros).sum(axis=1) == (mask | True).sum(axis=1)) | \
        ((mask.sum(axis=1) == 1) & row_suppress)

    
    col_check = (mask.sum(axis=0) == 0) | (mask.sum(axis=0) >= 2) | \
        ((mask | zeros).sum(axis=0) == (mask | True).sum(axis=0)) | \
        ((mask.sum(axis=0) == 1) & col_suppress) | \
        (axis >= 0)

    solved = row_check.all() & col_check.all() & suppressed.all().all()

    return solved, row_check, col_check


def _list_below(below):
    """
    Converts the list of cells below the value to a long form
    """
    need_supp = below.copy()
    
    # Sets up the index and column names ot make life easy
    if need_supp.index.names == [None]:
        need_supp.index.set_names('index', inplace=True)
    idx_names = need_supp.index.names

    if need_supp.columns.names == [None]:
        need_supp.columns.set_names('columns', inplace=True)
    col_names = need_supp.columns.names

    # Get a long form version of the table
    need_supp2 = need_supp.reset_index()
    need_supp2 = need_supp2.melt(id_vars=idx_names, )
    need_supp2 = need_supp2[need_supp2['value']].copy()

    # Casts the row and column names to a single index column with
    # tuples because otherwise this makes later work easier
    if (idx_names != ['index']) & (len(idx_names) == 1):
        need_supp2.rename(columns={idx_names[0]: 'index'}, inplace=True)
    elif (idx_names != ['index']):
        need_supp2['index'] = need_supp2[idx_names].apply(
            lambda x: tuple(x.values), axis=1)
    
    if (col_names != ['columns']) & (len(col_names) == 1):
        need_supp2.rename(columns={col_names[0]: 'columns'}, inplace=True)
    elif col_names != ['columns']:
        need_supp2['columns'] = need_supp2[col_names].apply(
            lambda x: tuple(x.values), axis=1)

    return need_supp2[['index', 'columns']].reset_index(drop=True)


def _make_coord_mask(index, column, below, row_optimum, col_optimum, **kwargs):
    """
    Builds all possible masks for hte data, optimizing the suppression
    """
    # Pulls out the columns that need to be solved. We start by checking
    # if the column value can be satisified by a matched "below" value
    # and if so, we use that as our default masking because it will
    # minimize the number of cells/masks we need to check. Otherwise, we
    # grab anything that hasn't been optimized.
    col_bel = below.loc[index].drop([column])
    if col_bel.any():
        col_sol = [col_bel.index[col_bel].values[0]]
    else:
        col_sol = (row_optimum).loc[index].drop([column])
        col_sol = col_sol.index[col_sol].values
    col_sol = [(index, cid) for cid in col_sol]
    
    # We repeat for the rows. If there's already a cell below, we pull
    # that value. Otherwise, we pull the values that need to be optimized
    row_bel = below[column].drop([index])
    if row_bel.any():
        row_sol= [row_bel.index[row_bel].values[0]]
    else:
        row_sol = (col_optimum)[column].drop([index])
        row_sol = row_sol.index[row_sol].values
    row_sol = [(rid, column) for rid in row_sol]
    
    # Builds the masks based on the supplied coordinates
    def _define_mask(coords):
        mask = below & False
        for (r, c) in coords:
            mask.loc[r, c] = True
        mask.loc[index, column] = True
        return mask

    # Builds the list of potential masks that solve the data
    if len(row_sol) == 0:
        choice_masks = [_define_mask([c]) for c in col_sol]
    elif len(col_sol) == 0:
        choice_masks = [_define_mask([c]) for c in row_sol]
    else:
        choice_masks = [_define_mask(c) 
                        for c in it.product(row_sol, col_sol)]

    return choice_masks


def _optimize_row_suppression_mask(row_below, col_below, squish_mask, squish_order, row_sum_supp, col_sum_supp, col_nz, penalty_mask, canidates, axis):
    """
    Optimizes the rows and columns that need to be suppressed

    We're looking for the optimum cell within a given row (or column if we
    transpose the data. We rank optimum rows in the following ways:
        1. Value is an avaliable squish value, based on the squish rank
        2. Value is being row suppressed
        3. Value is being column suppressed
        4. Column needs to be suppressed anwyay
        5. Column value is not zero
        6. Value is avaliable for suppression
    
    """
    # Optimizes the masking based on the columns contained in the row
    if axis >= 0:
        squish_order = squish_order.to_frame().values
    row_suppress_select = \
        (squish_mask) * (squish_order + 5) + \
        col_sum_supp * 4 + \
        col_below * 3 + \
        col_nz * 2 + \
        canidates * 1

    # Penalizes the penalty row where the penalty is half of the maximum
    penalty_suppressed = row_suppress_select - \
        (row_suppress_select.max(axis=1).to_frame().values + 1) / 2
    penalty_suppressed.mask((penalty_suppressed > row_suppress_select),
                            row_suppress_select,
                            inplace=True)
    penalty_suppressed.mask(penalty_suppressed < 0, 0, inplace=True)
    row_suppress_select.mask(penalty_mask, penalty_suppressed, 
                            inplace=True)

    avaliable = (row_below | row_sum_supp & canidates)
    row_suppress_select = (row_suppress_select * avaliable)
    row_best = row_suppress_select.max(axis=1).to_frame().values
    row_select = (row_suppress_select == row_best) & avaliable

    return row_suppress_select, row_select

