# Form definitions for user admin with flask-admin

from flask_security import current_user
from flask_admin.contrib.sqla import ModelView
from app import admin_obj, db
from security.userdb import User, Role


class AdminView(ModelView):
    def is_accessible(self):
        return current_user.has_role('Admin')

class UserView(AdminView):
    column_exclude_list = ('password')


admin_obj.add_view(UserView(User, db.session))
admin_obj.add_view(UserView(Role, db.session))
