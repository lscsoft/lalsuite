<!-- single time -->
<input type="radio" name="starttimevsgps" value="gps">GPS: <input type="text" size=10 maxlength=10 name="startgps">
            &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; or            &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
<input type="radio" checked="checked" name="starttimevsgps" value="time">Time: <select name="startmonth"><?php require '/usr1/ldbd/glue/var/php/script/form_month_list.php' ?></select>
                     <select name="startday"><?php require '/usr1/ldbd/glue/var/php/script/form_day_list.php' ?></select>
                     <select name="startyear"><?php require '/usr1/ldbd/glue/var/php/script/form_year_list.php' ?></select>
                     &nbsp&nbsp&nbsp
            <input type="text" name="starthour" value=00 size=2 maxlength=2>
            :<input type="text" name="startmin" value=00 size=2 maxlength=2>
            :<input type="text" name="startsec" value=00 size=2 maxlength=2>
                    &nbsp&nbsp&nbsp
            <select name="starttype">
                    <option value="24Hour">24 Hour</option>
                    <option value="am">am</option>
                    <option value="pm">pm</option>
            </select>
            <select name="startzone">
                    <option value="UTC">UTC</option>
                    <option value="Central">Central</option>
                    <option value="Pacific">Pacific</option>
            </select>
