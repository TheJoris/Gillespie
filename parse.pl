#!/usr/bin/perl

# perl script to parse schematic reaction files.

# Initialise all variables.  The large block of text below is 
# documentation printed when requested.

sub InitAll 
{
  $codename = "0";
  $reaction_file = "$codename.reactions";
  $component_file = "$codename.components";
  @name_list = ();
  @init_vals = ();
  @calc_vals = ();
  %components = ();
  $num_equi_blocks = -1;
  $num_prod_blocks = -1;
  $num_steps = 100000;
  $total_time_equi = 0;
  $total_time_run = 100;
  $time_step = 0.1;
  $volume = 1.0;
  $doubling_time = -1.; # no growth at all
#   $doubling_time_std = 0; # exact doubling time
  $help = <<EOF;

Parses schematic reaction files and generates an input file for Gillespie
on stdout.  Schematic reaction files are parsed line by line, each line is
one of the following forms.

1. Blank lines and lines not containing ":" are ignored.

2. All text after '#' is ignored.

3. A line of the type
 <value> : <identifier> 
resets default values.  The <identifier> is from the following table

-------------------------------------------------------------------------
 <identifier>          minimum components in <identifier>  default
-------------------------------------------------------------------------
 number of equilibrium block   "num" "eq"                   -1
 number of production blocks   "num" "prod"                 -1
 number of steps               "num" "steps"                100000
 maximum total time for equi.  "total" "eq"                 0
 maximum total time for prod.  "total" "prod"               100
 time step of write out        "time"                       0.1 
 code name                     "name"                       0
 reaction file                 "react" "file"               0.reactions
 component file                "comp" "file"                0.components
 volume factor                 "vol"                        1
 doubling time                 "doubling"                   -1
 doubling time noise           "std"                        0
 global degradation rate       "degradation"                0.
-------------------------------------------------------------------------

The simulation runs until either the number of simulation steps (this is "number of
blocks" times "number of steps") or the total time has been reached. Both values can
be set independently for the equilibrium and the production phase of the algorithm.

Changing the code name also has the effect of changing the file names.
Setting the volume factor only affects reaction constants and initial conditions 
below the respective line, so it should be one of the first things to set.
A negative doubling time defaults to no growth.

4. A line of the type
 <space_separated_list_of_components> : components
declares components, and is used to determine a particular order of
appearance, otherwise components appear in the order in which they
appear in reactions.

5. A line of the type
 <space_separated_list_of_components> : constant
also declares components, which are now to be taken as constant during
cell division.

6. A line of the type
 <component> <equals> <val> : initial 
causes <component> name to have a initial value <val> (default 0).  
The separator <equals> can be anything.

7. A line of the type
 <component> = <comp1> + <comp2> [ + <compn> ] : calc
causes <component> to equal the sum of all the other components after each
Gillespie step.

8. A line of the type 
 <reactants> <transform> <products> : <k_values>
declares a reaction with the following grammar.

Production reactions (can write "null" or "NULL"),
 null --> A  : k = <val>
 null --> A + B : k = <val>

Annihilation reactions (can write "null" or "NULL"),
 A --> null : k = <val>
 A + B --> null : k = <val>

Simple reactions,
 A + B --> C + D + E : k = <val>

Production reactions, in which the reactants are not destroyed,
 A + B +--> C + D : k = <val>
is equivalent to
 A + B --> A + B + C + D : k = <val>

Equilibriums, which are converted to reaction pairs,
 A + B <--> C + D : k = <val_forward> : k = <val_backward>
is equivalent to
 A + B --> C + D : k = <val_forward>
 C + D --> A + B : k = <val_backward>

Finally simple and production reactions can be written with reactants
and products as lists of alternatives.  In such cases all possible
pairs of alternatives are implemented, with the same k, thus:
 A | B | C --> null : k = <val_1>
 A + B +--> C | D : k = <val_2>
is equivalent to
 A --> null : k = <val_1>
 B --> null : k = <val_1>
 C --> null : k = <val_1>
 A + B --> A + B + C : k = <val_2>
 A + B --> A + B + D : k = <val_2>

9. Instead of using the single arrow, a double arrow denotes a delayed reaction:
 A + B ==> C + D : k = <val_k> : T = <val_T> : s = <val_s>
The delays will be Gaussian distributed with mean <val_T> and standard deviation <val_s>.

10. Instead of a constant value for the reaction constant k, the following
may be used to denote Hill function dependence:
 A --> B : k = Hill(k0,X,K,n)
which translates to
k = k0 * [X]^n / ( [X]^n + K^n )
Please note that repression can be modelled using negative exponents n.

At the end, reaction and component files suitable for direct reading
by Gillespie code are generated.  Standard output can be redirected to
Gillespie.inp (the default avoids accidental overwriting).  Various
comments are generated on standard error.

EOF
}


# Returns id number of a component of a given name, either one
# previously used or a new one which is initialised appropriately.
sub component_id 
{
  local ( $name, $id );
  $name = $_[0];
  $name =~ s/\s//g;

  # check, if component exists
  if( exists $components{$name} )
  {
  	$id = $components{$name};
  } 
  else # add new component
  {
    $id = scalar(@name_list);
    $components{$name} = $id;
    @init_vals[$id] = 0;
    @comp_const[$id] = 0;
    @name_list[$id] = $name;
  }
  return $id;
}


# Add a list of components, in the order in which they appear.
sub add_components 
{
  local ( $name );
  @names = split(' ',$_[0]);
  foreach $name (@names)
  {
    component_id($name);
  }
}


# Add a list of components, in the order in which they appear and sets them to constant
sub add_components_const
{
  local ( $name );
  @names = split(' ',$_[0]);
  foreach $name (@names)
  {
    $id = component_id($name);
    @comp_const[$id] = 1
  }
}


# Set or reset the initial concentration for a given component.
sub set_initial 
{
  local ( $name, $expr, $val, $id );
  ( $name, $expr, $val ) = split(' ',$_[0]);
  
  $id = component_id($name);

  # check, if component should be constant or if it is effected by volume
  if( not @comp_const[$id] )
  {
    $val = $val * $volume;
  }
  printf STDERR "Initialising component %s to %i\n", $name, $val;
  @init_vals[$id] = $val;
}


# Set or reset the calculation for the concentration for a given component.
sub set_calculation
{
  local ( $name, $comps_str, $component, $i, @id_list, @comps );
  ( $name, $comps_str ) = split('[=]+',$_[0]);
  
  $id = component_id($name);

  # remember all the components going into the sum
  @comps = split('\+', $comps_str);
  @id_list = ();
  $i = 0;
  foreach $component (@comps)
  {
    $c_id = component_id($component);
    $id_list[$i] = $c_id;
    $i = $i + 1;
  }
  
  # set the list for this component
  @comp_const[$id] = 2;
  @calc_vals[$id] = [@id_list];
  
  printf STDERR "Component %s is calculated from %s\n", $name, $comps_str;
}


# Add a reaction, specifically compute the output lines for later
# writing.
sub add_reaction 
{
  local ( $reacts, $prods, $kval, $component, $id, $nreacts, $nprods );
  ( $reacts, $prods, $kval, $time, $sigma ) = ( $_[0], $_[1], $_[2], $_[3], $_[4] );

  # check kval for special identifiers

  if ($kval =~ /Hill\((.*)\)/)
  {
    # found Hill function type ==> extract values
    ( $kval, $HillComp, $HillConst, $HillCoeff ) =  split(",",$1,4);
    $HillComp = component_id( $HillComp );
    $HillConst =~ s/\s//g;
    $HillCoeff =~ s/\s//g;
    $kstr = sprintf( "Hill(%f,%s,%f,%f)", $kval, @name_list[$HillComp], $HillConst, $HillCoeff );
  }
  else # no special identifier ==> strip
  {
    ( $HillComp, $HillConst, $HillCoeff ) = (0,0,0);
    $kstr = sprintf( "%f", $kval ); #TODO: change such that arithmic expression is conv. to number.
  }

  if( $time == 0 && $sigma == 0)
  {
    printf STDERR "Adding reaction %s --> %s, k = %s\n", $_[0], $_[1], $kstr;
  }
  else
  {
    printf STDERR "Adding reaction %s ==> %s, k = %s, T = %f, s = %f\n", $_[0], $_[1], $kstr, $_[3], $_[4];
  }



  # strip constants
  $kval =~ s/\s//g;
  $time =~ s/\s//g;
  $sigma =~ s/\s//g;

  # init output line
  $line2 = "";

  # look for reactants
  if ($reacts =~ /null/ || $reacts =~ /NULL/)
  {
    @reacts = ();
    $nreats = 0;
    $line2 .= "0 ";
  }
  else
  {
    @reacts = split(/\+/,$reacts);
    $nreacts = scalar(@reacts);
    $addmore = 0;
    foreach $component (@reacts)
    {
      if( $addmore )
      {
        $line2 .= "+ ";
      }
      $id = component_id($component);
      $bit = sprintf("X %i ",$id);
      $line2 .= $bit;
      $addmore = 1;
    }
  }

  $line2 .= "->";

  # add products
  if ($prods =~ /null/ || $prods =~ /NULL/)
  {
    $line2 .= " 0";
    $nprods = 0;
  }
  else
  {
    @prods = split(/\+/,$prods);
    $nprods = scalar(@prods);
    $addmore = 0;
    foreach $component (@prods)
    {
      if ($addmore)
      {
        $line2 .= " +";
      }
      $component =~ s/\s//g;
      @data = split(/\*/,$component);
      if( scalar(@data) == 1 )
      {
        $id = component_id($component);
        $line2 .= sprintf(" 1 X %i",$id);
      }
      elsif( scalar(@data) == 2 )
      {
        $id = component_id(@data[1]);
        $line2 .= sprintf(" %s X %i",@data[0],$id);
      }
      else
      {
        print "Ill-formed reaction line.";
        exit(0);
      }
      $addmore = 1;
    }
  }

  # possibly modify the $kval according to current volume factor
  $HillConst = $HillConst * $volume;
  if( 0 == $nreacts )
  {
    $kval = $kval * $volume;
  }
  if( 2 == $nreacts )
  {
    $kval = $kval / $volume;
  }

  # set the parameter output line
  $line1 = sprintf(
            "%f\t%i\t%i\t%f\t%f\t%i\t%f\t%f\tRateConstantk_Nreactants_Nproducts_Time_Sigma_HillComp_HillConst_HillCoeff",
            $kval,$nreacts,$nprods,$time,$sigma,$HillComp,$HillConst,$HillCoeff
           );

  # save the lines
  push(@react1lines, $line1);
  push(@react2lines, $line2);
}


# iterates over all species and adds a degradation reaction
sub add_degradation_reactions
{
  local ( $degradation, $ncomp, $id );
  $degradation = $_[0];

  $ncomp = scalar(@name_list);
  for ($id=0; $id<$ncomp; $id++)
  {
    if( 0 == $comp_const[$id] )
    {
      add_reaction( @name_list[$id], "null", $degradation );
    }
  }
}


# Parse a single line from a schematic reaction file.
sub parseline 
{
  local ( $line, $front, $back, $reacts, $prods );
  $line = $_[0];
  # strip comments
  $line =~ s/#.*//;
  
  if ($line =~ /:/)
  {
    ( $front, $back ) = split(":",$line,2);
    # look for variables
    if ($back =~ /name/)
    {
      $codename = $front;
      $codename =~ s/\s//g;
      $reaction_file = "$codename.reactions";
      $component_file = "$codename.components";
      print STDERR "Setting codename to $codename\n";
    }
    elsif ($back =~ /vol/)
    {
      $volume = $front;
      $volume =~ s/\s//g;
      print STDERR "Setting volume factor to $volume\n";
    }
    elsif ($back =~ /file/)
    {
      if ($back =~ /react/)
      {
        $reaction_file = $front; $reaction_file =~ s/\s//g;
        print STDERR "Setting reaction file to $reaction_file\n";
      }
      elsif ($back =~ /comp/)
      {
        $component_file = $front; $component_file =~ s/\s//g;
        print STDERR "Setting component file to $component_file\n";
      }
    }
    elsif ($back =~ /num/)
    {
      if ($back =~ /eq/)
      {
        $num_equi_blocks = $front;
        $num_equi_blocks =~ s/\s//g;
        printf STDERR "Setting # equilibration blocks to %i\n", $num_equi_blocks;
      }
      elsif ($back =~ /prod/)
      {
        $num_prod_blocks = $front;
        $num_prod_blocks =~ s/\s//g;
        printf STDERR "Setting # production blocks to %i\n", $num_prod_blocks;
      }
      elsif ($back =~ /step/)
      {
        $num_steps = $front;
        $num_steps =~ s/\s//g;
        printf STDERR "Setting # steps to %i\n", $num_steps;
      }
    }
    elsif ($back =~ /doubling/)
    {
      $doubling_time = $front;
      $doubling_time =~ s/\s//g;
      printf STDERR "Setting the doubling_time to %f\n", $doubling_time;
    }
    elsif ($back =~ /std/)
    {
      $doubling_time_std = $front;
      $doubling_time_std =~ s/\s//g;
      printf STDERR "Setting the doubling_time_std to %f\n", $doubling_time_std;
    }
    elsif ($back =~ /time/)
    {
      $time_step = $front;
      $time_step =~ s/\s//g;
      printf STDERR "Setting the time_step to %f\n", $time_step;
    }
    elsif ($back =~ /total/)
    {
      if ($back =~ /eq/)
      {
        $total_time_equi = $front;
        $total_time_equi =~ s/\s//g;
        printf STDERR "Setting total time of the equilibrium run to %f\n", $total_time_equi;
      }
      elsif ($back =~ /run/ or $back =~ /prod/)
      {
        $total_time_run = $front;
        $total_time_run =~ s/\s//g;
        printf STDERR "Setting total time of the production run to %f\n", $total_time_run;
      }
      else
      {
        $total_time_run = $front;
        $total_time_run =~ s/\s//g;
        $total_time_equi = $total_time_run;
        printf STDERR "Setting total time of the equilibrium and the production run to %f\n", $total_time_run;
      }
    }
    elsif ($back =~ /init/)
    {
      set_initial($front);
    }
    elsif ($back =~ /calc/)
    {
      set_calculation($front);
    }
    elsif ($back =~ /comp/)
    {
      add_components($front);
    }
    elsif ($back =~ /constant/)
    {
      add_components_const($front);
    }
    elsif ($back =~ /degradation/)
    {
      $degradation = $front;
      $degradation =~ s/\s//g;
      printf STDERR "Setting degradation rate for all species to %f\n", $degradation;
      add_degradation_reactions( $degradation );
    }
    # look for reactions
    elsif ($front =~ /<[=-]+>/)   # <--> or <==>
    {
      ( $optreacts, $optprods ) = split('<[=-]+>',$front);
      @optreacts = split(/\|/,$optreacts);
      @optprods = split(/\|/,$optprods);
      if ($back =~ /:/)
      {
        ( $kforward, $kback, $time, $sigma ) = split(/:/,$back,4);
      }
      else
      {
        $kforward = $back;
        $kback = $back;
        $time = 0.;
        $sigma = 0.;
      }
      $kforward =~ s/.*=\s*//;
      $kback =~ s/.*=\s*//;
      $time =~ s/.*=\s*//;
      $sigma =~ s/.*=\s*//;
      foreach $reacts (@optreacts)
      {
        foreach $prods (@optprods)
        {
          add_reaction($reacts,$prods,$kforward,$time,$sigma);
          add_reaction($prods,$reacts,$kback,$time,$sigma);
        }
      }
    }
    elsif ($front =~ /\+[=-]+>/)  # +--> or +==>
    {
      ( $optreacts, $optprods ) = split('\+[=-]+>',$front);
      @optreacts = split(/\|/,$optreacts);
      @optprods = split(/\|/,$optprods);
      if ($back =~ /:/)
      {
        ( $kval, $time, $sigma ) = split(/:/,$back,3);
      }
      else
      {
        $kval = $back;
        $time = 0.;
        $sigma = 0.;
      }
      $kval =~ s/.*=\s*//;
      $time =~ s/.*=\s*//;
      $sigma =~ s/.*=\s*//;
      foreach $reacts (@optreacts)
      {
        foreach $prods (@optprods)
        {
          add_reaction($reacts,"$reacts + $prods",$kval,$time,$sigma);
        }
      }
    }
    elsif ($front =~ /[=-]+>/)  # --> or ==>
    {
      ( $optreacts, $optprods ) = split('[=-]+>',$front);
      @optreacts = split(/\|/,$optreacts);
      @optprods = split(/\|/,$optprods);
      if ($back =~ /:/)
      {
        ( $kval, $time, $sigma ) = split(/:/,$back,3);
      }
      else
      {
        $kval = $back;
        $time = 0.;
        $sigma = 0.;
      }
      $kval =~ s/.*=\s*//;
      $time =~ s/.*=\s*//;
      $sigma =~ s/.*=\s*//;
      foreach $reacts (@optreacts)
      {
        foreach $prods (@optprods)
        {
          add_reaction($reacts,$prods,$kval,$time,$sigma);
        }
      }
    }
    else
    {
      printf STDERR "Didn't recognise %s\n",$line;
    }
  }
}


# Parse a schematic reaction file line by line.
sub ParseFile 
{
  $in_file = $_[0];
  open(FIN, $in_file) || die "Cannot open $in_file: $!";
  printf STDERR "Reading from $in_file\n";
  while ($line = <FIN>)
  {
    chop($line);
    parseline($line);
  }
  close(FIN);
}


# Write a list of current reactions to a file using PRs format.
sub WriteReactionFile 
{
  $out_file = $_[0];
  open(FOUT, "> $out_file") || die "Cannot open $out_file: $!";
  $nreact = scalar(@react1lines);
  printf FOUT "%i\t\tNumber_of_reaction_channels\n\n",$nreact;
  for ($i=0; $i<$nreact; $i++)
  {
  	printf FOUT "%s\n%s\n",@react1lines[$i],@react2lines[$i];
  }
  print FOUT "\n";
  close(FOUT);
  printf STDERR "Written %i reactions to $out_file\n",$nreact;
}


# Write a list of current components to a file using PRs format.
sub WriteComponentFile 
{
  local ( $ncomp, $id, $c_id, @c_ids );
  $out_file = $_[0];
  open(FOUT, "> $out_file") || die "Cannot open $out_file: $!";
  $ncomp = scalar(@name_list);
  printf FOUT "%i\n",$ncomp;
  for ($id=0; $id<$ncomp; $id++)
  {
    printf FOUT "%i\t\t%s\t%d", @init_vals[$id], @name_list[$id], @comp_const[$id];
    if( @comp_const[$id] == 2 )
    {
      @c_ids = @{ $calc_vals[$id]};
      printf FOUT "\t%i", scalar (@c_ids);
      foreach $c_id (@c_ids)
      {
        printf FOUT "\t%i", $c_id;
      }
    }
    printf FOUT "\n"
  }
  print FOUT "\n";
  close(FOUT);
  printf STDERR "Written %i components to $out_file\n",$ncomp;
}


# Write a control file to stdout using current information, suitable
# for PRs Gillespie code.  It is intended this should be re-directed
# as appropriate.
sub WriteGillespieInp 
{
  printf "%s\t\t\tName\n",$codename;
  printf "%i\t\t\tNumber_of_components\n",scalar(@name_list);
  printf "%i\t\t\tNumber_of_reactions\n",scalar(@react1lines);
  printf "%i\t%i\t%i\tNumber_eq_blocks;_Number_prod_blocks;_Number_steps\n",
          $num_equi_blocks,$num_prod_blocks,$num_steps;
  printf "%f\t\t\tTime_step\n",$time_step;
  printf "%f\t%f\t\tTotal_time\n\n",$total_time_equi,$total_time_run;
  printf "%f\t%f\t\tDoubling_time\n\n", $doubling_time, $doubling_time_std;
}

# This is where the action starts.

InitAll;

if( scalar(@ARGV) == 0 )
{
  printf STDERR "Usage: ./parse.pl <input_files> > Gillespie.inp\n";
  printf STDERR "./parse.pl --help for more info\n";
  exit(0);
}

if( @ARGV[0] =~ /--help/ )
{
  print "Usage: ./parse.pl <input_files> > Gillespie.inp\n";
  print $help;
  exit(0);
}

foreach $file (@ARGV) 
{
  ParseFile($file);
}

WriteReactionFile($reaction_file);
WriteComponentFile($component_file);
WriteGillespieInp;

# End of file.
