# The \`compartmentalFit\` object

The object returned by \[fitCompartmentalModel()\]. Its components are
public, but the accessors (see \[compfit-accessors\]) provide a stable,
validated interface and are the recommended way to reach them.

## Components

- \`method\`, \`success\`, \`error_msg\`:

  Run metadata.

- \`point\`:

  \*(point methods only)\* \`list(initial_state, parms)\` on the natural
  scale.

- \`best_state\`:

  Environment with \`\$error\` (best loss value reached) and \`\$par\`
  (best normalised parameter vector).

- \`fit_raw\`:

  Raw optimiser/sampler return (\`optim\`, \`DEoptim\`, or a Turing
  chain list), or a \`fitError\` object on failure.

- \`samples\`:

  \*(bayes only)\* list with: \`draws\` (one row per joint MCMC sample;
  fitted parameters as normalised \`name_n\` columns, hyperparameters on
  natural scale — use \[posterior_draws()\] to convert all to natural
  scale), \`summary\` (\`parameter, mean, sd, rhat, ess\`),
  \`prior_spec\`, \`like_specs\`, \`model_code\` (generated Turing
  \`@model\` source), \`scale_cols\`, \`stream_names\`, \`n_chains\`,
  \`iter\`.

- \`loss\`:

  The loss closure (dropped on \`save_fit()\`, rebuilt by \[load_fit()\]
  on load).

- \`model\`:

  List: \`structure\`, \`sap\`, \`expressions\`,
  \`compartmental_function\` (R closure), \`julia_code\` (character),
  \`date\`, and \`modelParams\` (the source sheet).

- \`data\`:

  Prepared data; all matrices are years x streams: \`data_combined\`,
  \`matrix_data_points\`, the observation/censoring masks (\`obs_mask\`;
  \`cens_mask\`/\`limit_mat\` for left-censored \`\<L\`/\`\<=L\`;
  \`lcens_mask\`/\`llimit_mat\` for right-censored \`\>L\`/\`\>=L\`),
  \`weight_matrix\`, \`average_matrix\`, \`names_data_points\`, and
  \`cumulative_cols\` (stream indices differenced from cumulative to
  annual).

- \`time_grid\`:

  List: \`time\`, \`startpoint\`, \`endpoint\`, \`partition\`,
  \`cutoff\`.

- \`bounds\`:

  List: \`lower\`, \`upper\`, \`init_norm\`.

- \`solver\`:

  The \[solver_control()\] list used.

- \`meta\`:

  List: \`checkpoint_file\`.

## See also

\[fitCompartmentalModel()\], \[compfit-accessors\], \[save_fit()\],
\[load_fit()\], \[summary.compartmentalFit()\],
\[print.compartmentalFit()\]
