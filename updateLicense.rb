#!/usr/bin/env ruby
# This script adds/update the license header of bogus source files
# Any copyright is dedicated to the Public Domain.
# http://creativecommons.org/publicdomain/zero/1.0/ 

GIT_CMD =  'git ls-files "*.cpp" "*.hpp" "*.h" '
EXCLUDE = []

COPYRIGHT_LINES = [
  'Copyright 2016 Gilles Daviet <gdaviet@gmail.com>'
]

MPL_TEXT = <<eol
 * This file is part of bogus, a C++ sparse block matrix library.
 *
 * {COPYRIGHT}
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
eol

GPL_TEXT = <<eol
 * This file is part of So-bogus, a C++ sparse block matrix library and
 * Second Order Cone solver.
 *
 * {COPYRIGHT}
 *
 * So-bogus is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.

 * So-bogus is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with So-bogus.  If not, see <http://www.gnu.org/licenses/>.
eol

CC0_TEXT = <<eol
 * Any copyright is dedicated to the Public Domain.
 * http://creativecommons.org/publicdomain/zero/1.0/ 
eol

class Licenser
  def initialize( text )
    expand_license!( text )
  end

  def expand_license!( text )
    @text = ''
    text.each_line do |l|
      if l.include? '{COPYRIGHT}' then
        COPYRIGHT_LINES.each { |c|
          @text << l.sub( '{COPYRIGHT}', c)
        }
      else
        @text << l
      end
    end
  end
  
  def add_copyright out
    out << "/*" 
    out << @text
    out << "*/" 
  end

  def process path
    puts "Processing #{path}\n"

    out_lines = []

    comment_start = '/*'
    comment_end = '*/'

    inside_comment = false
    done = false
    
    File.open( path ) do |f|
      f.each_line do |l|

        line = l.strip
            
        inside_comment = l.start_with?( comment_start ) unless done or inside_comment

        if inside_comment
          done = line.end_with?( comment_end ) 
          inside_comment = ! done
          add_copyright( out_lines ) if done
        else 

          unless done then

            done = ! ( inside_comment or line.empty? )
            if done then
              add_copyright( out_lines )
              out_lines << "" # add an extra empty line
            end
          end

          out_lines << l unless inside_comment 

        end


      end
    end

    File.open( path, "w" ) do |f|
      out_lines.each do |l|
        f.puts l
      end
    end

  end

end


license = ARGV[0].downcase if ARGV[0]
case license
  when "gpl"
    Text = GPL_TEXT
  when "mpl"
    Text = MPL_TEXT
  when "cc0"
    Text = CC0_TEXT
  else
    raise 'Unknown license'
end

licenser = Licenser.new( Text )

(1..(ARGV.count-1)).each do |i| 
  prefix = ARGV[i] 
  IO.popen(GIT_CMD) do |f|
    f.each_line do |l|
      l.strip!
      ok = l.start_with?( prefix ) 
      EXCLUDE.each{ |str| ok &&= !l.include?( str ) } if ok
      if ok then
        licenser.process l
      end
    end
  end
end


